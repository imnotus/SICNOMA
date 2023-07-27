#include <bits/stdc++.h>
//#include "MT.h"

using namespace std;
float inf = numeric_limits<float>::infinity();
using ll = long long;
const double radius = 100;
const double GAMMA = pow(10, 6.0 / 10.0);
const double end_time = 100000;
constexpr double PI = 3.14159265358979323846264338;
const double theta = pow(10, 0.4);
//const double delta = 2.0 / pass_loss_exponent;
const double lambda_BS = 1.;
//const double POWER = 1.0;

//出力ファイルを開く
ofstream outputfile1;
ofstream outputfile2;
ofstream outputfile;

struct pos {
    double x;
    double y;
};

struct terminal {
    pair<double, double> pos;
    //double to_BS[num_BS];
    double dst_to_origin;
    int nearest_BS;
    double level;
    double coef;
    bool state; //True: in , False: out
    bool operator<( const terminal& right ) const {
            return coef == right.coef ? true : coef < right.coef;
    }
};


double urand(){
    double m, a;
    m = RAND_MAX + 1.0;
    a = (rand() + 0.5)/m;
    a = (rand() + a)/m;
    return (rand() + a)/m;
}

//乱数発生
double my_rand(double Min, double Max) {
    mt19937 mt{ random_device{}() };
    //random_device rd;
    //default_random_engine eng(rd());
    uniform_real_distribution<double> distr(Min, Max);
    return distr(mt);
}

//円形領域内の座標をランダムに取得
pair<double, double> coordinate() {
    //init_genrand(i);
    double r = radius * sqrt(urand()); //[0, 1)
    double theta = PI * (2. * urand() - 1.); //(-1, 1)
    while (theta < -PI || PI < theta) {
        theta = PI * (2. * urand() - 1.);
    }
    double x = r * cos(theta);
    double y = r * sin(theta);
    pair<double, double> pos = make_pair(x, y);
    return pos;
}

//2点間の距離を計算
double cal_dst(pair<double, double> pos1, pair<double, double> pos2) {
    return sqrt((pos1.first - pos2.first)*(pos1.first - pos2.first) + (pos1.second - pos2.second)*(pos1.second - pos2.second));
}

double distance(pos pos1, pos pos2) {
    return sqrt((pos1.x - pos2.x) * (pos1.x - pos2.x) + (pos1.y - pos2.y) * (pos1.y - pos2.y));
}

//正規分布乱数発生
double normrand(){
    return sqrt(-2*log(urand()))*cos(2*M_PI*urand());
}


double gauss_rand(double mu, double sig) {
    //mt19937 mt{ std::random_device{}() };
    random_device rd;
    default_random_engine eng(rd());
    normal_distribution<> dist(mu, sig);
    double g = dist(eng);
    return g;
}


double exp_dist(double lambda) {
    //double g = genrand_real3();
    double g = urand();
    double tau = - log(1 - g) / lambda;
    return tau;
}

double poisson_dist(double lambda) {
    double sum = 0, k;
    for (k = 0; sum < 1; k++) {
        sum += exp_dist(lambda);
    }
    return  k;
}

//ポアソン乱数
double poisson(double lambda) {
    random_device rd;
    default_random_engine eng(rd());
    poisson_distribution<> dist(lambda);
    return dist(eng);
}

//NOMA
void SIC_NOMA_thp(double lambda_IoT, double alpha, double noise) {
    double success = 0.0;
    cout << "pass loss exp : " << alpha << endl;
    
    //座標読み込み
    string filename1 = "BS_pos.txt";
    vector<pair<double, double>> BS_pos;
    ifstream readingfile;
    readingfile.open(filename1);
    ifstream ifs("BS_pos.txt");
        if (ifs.fail()) {
           cerr << "Cannot open file\n";
           exit(0);
        }
    double x, y;
    while (ifs >> x >> y) {
        BS_pos.push_back(make_pair(x, y));
        //cout << x << " " << y << endl;
    }
    ifs.close();
    unsigned long num_BS = BS_pos.size();
//    int num_BS = poisson_dist(PI * radius * radius * lambda_BS);
//    //BSの座標を設定
//    vector<pair<double, double>> BS_pos(num_BS);
    pair<double, double> origin = make_pair(0, 0);
    BS_pos.at(0) = origin;
//    for (int i = 1; i < num_BS; i++) {
//        BS_pos.at(i) = coordinate();
//        outputfile1 << BS_pos.at(i).first << " " << BS_pos.at(i).second << endl;
//    }
    
    
    bool p_flag[4] = {true, true, true, true};
    cout << "Progress is 0%";
    for (int t = 0; t < end_time; t++) {
        int num_IoT = poisson_dist(PI * radius * radius * lambda_IoT);
        
        double SI = 0;
        vector<terminal> device(num_IoT);
        vector<int> acl;
        //端末情報を初期化
        for (int i = 0; i < num_IoT; i++) {
            //pair<double, double> pos = coordinate();
            device.at(i).pos = coordinate();
            double dst_to = cal_dst(device.at(i).pos, origin);
            device.at(i).dst_to_origin = dst_to;
            if (dst_to < 30) { //原点基地局にアクセスする端末を探す
                for (int j = 0; j < num_BS; j++) {
                    double dst = cal_dst(device.at(i).pos, BS_pos.at(j));
                    device.at(i).state = true; //送信先が原点ならtrue
                    if (dst < device.at(i).dst_to_origin) {
                        device.at(i).state = false;
                        break;
                    }
                    if (j == num_BS - 1) acl.push_back(i);
                }
            }
            
            //端末-基地局間のフェージング係数を設定
            double H = exp_dist(1.0);
            double coef = H / pow(device.at(i).dst_to_origin, alpha);
            device.at(i).coef = coef;
            SI += coef;
        }
        
        pair<double, double> FD_pos;
        pair<double, double> SD_pos;
        int s = (int)acl.size();
        double big = 0, big2 = 0;
        for (int i = 0; i < s; i++) {
            double P = device.at(acl.at(i)).coef;
            if (big2 < P) {
                if (big < P) {
                    big = P;
                    FD_pos = device.at(acl.at(i)).pos;
                } else {
                    big2 = P;
                    SD_pos = device.at(acl.at(i)).pos;
                }
            }
        }
        
        SI -= big;
        double SINR = big / (SI + noise);
        if (SINR > theta) {
            success++;
            outputfile1 << FD_pos.first << " " << FD_pos.second << endl;
        }
        SI -= big2;
        SINR = big2 / (SI + noise);
        if (SINR > theta) {
            outputfile2 << SD_pos.first << " " << SD_pos.second << endl;
            //cout << IoT_pos.first << " " << IoT_pos.second << endl;
            success++;
        }
        
        
        //進捗状況を表示
        double progress = t / end_time;
        if (progress > 0.8 && p_flag[3]) {cout << "...80%"; p_flag[3] = false;}
        else if (progress > 0.6 && p_flag[2]) {cout << "...60%"; p_flag[2] = false;}
        else if (progress > 0.4 && p_flag[1]) {cout << "...40%"; p_flag[1] = false;}
        else if (progress > 0.2 && p_flag[0]) {cout << "...20%"; p_flag[0] = false;}
    }
    cout << "...100%" << endl << endl;
    
    //double res = success / end_time;
    //double I = thp_theo(lambda_IoT, lambda_BS, alpha);
    //cout << "Throughput : " <<  res << " " << I << endl << endl;
    //outputfile << res << " " << I << " ";
}


//Power allocation
void PA_NOMA_pos(double lambda_IoT, double alpha, double noise, double L) {
    double success = 0.0;
    cout << "pass loss exp : " << alpha << endl;
    
    //座標読み込み
    string filename0 = "BS_pos.txt";
    vector<pair<double, double>> BS_pos;
    ifstream readingfile;
    readingfile.open(filename0);
    ifstream ifs("BS_pos.txt");
        if (ifs.fail()) {
           cerr << "Cannot open file\n";
           exit(0);
        }
    double x, y;
    while (ifs >> x >> y) {
        BS_pos.push_back(make_pair(x, y));
        //cout << x << " " << y << endl;
    }
    ifs.close();
    unsigned long num_BS = BS_pos.size();
    cout << num_BS << endl;
//    int num_BS = poisson_dist(PI * radius * radius * lambda_BS);
//    //BSの座標を設定
//    vector<pair<double, double>> BS_pos(num_BS);
    pair<double, double> origin = make_pair(0, 0);
    BS_pos.at(0) = origin;
//    for (int i = 1; i < num_BS; i++) {
//        BS_pos.at(i) = coordinate();
//        outputfile1 << BS_pos.at(i).first << " " << BS_pos.at(i).second << endl;
//    }
    
    outputfile << "# L = " << L << endl;
    bool p_flag[4] = {true, true, true, true};
    cout << "Progress is 0%";
    for (int t = 0; t < end_time; t++) {
        int num_IoT = poisson_dist(PI * radius * radius * lambda_IoT);
        
        double SI = 0;
        vector<double> LEV(L, 0);
        vector<terminal> device(num_IoT);
        vector<int> acl;
        //端末情報を初期化
        for (int i = 0; i < num_IoT; i++) {
            //pair<double, double> pos = coordinate();
            device.at(i).pos = coordinate();
            double dst_to = cal_dst(device.at(i).pos, origin);
            device.at(i).dst_to_origin = dst_to;
            if (dst_to < 5) { //原点基地局にアクセスする端末を探す
                for (int j = 0; j < num_BS; j++) {
                    double dst = cal_dst(device.at(i).pos, BS_pos.at(j));
                    device.at(i).state = true; //送信先が原点ならtrue
                    if (dst < device.at(i).dst_to_origin) {
                        device.at(i).state = false;
                        break;
                    }
                    if (j == num_BS - 1) acl.push_back(i);
                }
            }
            
            if (device.at(i).state) { //原点にアクセスするならallocateされた電力を足す．
                int level = L * device.at(i).dst_to_origin;;
                if (level == 0) level = 1;
                else if (level >= L) level = L;
                device.at(i).level = level;
                device.at(i).coef = theta * pow(theta + 1, L - level);
                LEV.at(level-1) += theta * pow(theta + 1, L - level);
            } else {
                //端末-基地局間のフェージング係数を設定
                double H = exp_dist(1.0);//gauss_rand(0, 1);
                double coef = H / pow(device.at(i).dst_to_origin, alpha);
                device.at(i).coef = coef;
                SI += coef;
            }
        }
        
        int s = (int)acl.size();
        for (int i = 0; i < s; i++) {
            int f = acl.at(i);
            double P = device.at(f).coef;
            int l = device.at(f).level;
            double si = 0;
            for (int j = 0; j <= l; j++) {
                si += LEV.at(j);
            }
            double SINR = P / (SI + si - P + noise);
            if (SINR > theta) {
                success++;
                outputfile << device.at(i).pos.first << " " << device.at(i).pos.second << endl;
            }
        }
        
        
        //進捗状況を表示
        double progress = t / end_time;
        if (progress > 0.8 && p_flag[3]) {cout << "...80%"; p_flag[3] = false;}
        else if (progress > 0.6 && p_flag[2]) {cout << "...60%"; p_flag[2] = false;}
        else if (progress > 0.4 && p_flag[1]) {cout << "...40%"; p_flag[1] = false;}
        else if (progress > 0.2 && p_flag[0]) {cout << "...20%"; p_flag[0] = false;}
    }
    cout << "...100%" << endl << endl;
    outputfile << endl << endl;
    
    //double res = success / end_time;
    //double I = thp_theo(lambda_IoT, lambda_BS, alpha);
    //cout << "Throughput : " <<  res << " " << I << endl << endl;
    //outputfile << res << " " << I << " ";
}

//Power allocated NOMA
void PA_NOMA_thp(double lambda_IoT, double alpha, double noise, double L) {
    double success = 0.0;
    cout << "pass loss exp : " << alpha << endl;
    
    outputfile << lambda_IoT << " ";
    bool p_flag[4] = {true, true, true, true};
    cout << "Progress is 0%";
    for (int i = 0; i < end_time; i++) {
        int num_IoT = poisson_dist(PI * radius * radius * lambda_IoT);
        int num_BS = poisson_dist(PI * radius * radius * lambda_BS);
        //BSの座標を設定
        vector<pair<double, double>> BS_pos(num_BS);
        pair<double, double> origin = make_pair(0, 0);
        BS_pos.at(0) = origin;
        for (int i = 1; i < num_BS; i++) {
            BS_pos.at(i) = coordinate();
        }
        
        double SI = 0;
        vector<double> LEV(L, 0);
        vector<terminal> device(num_IoT);
        vector<int> acl;
        //端末情報を初期化
        for (int i = 0; i < num_IoT; i++) {
            pair<double, double> pos = coordinate();
            double dst_to = cal_dst(pos, origin);
            device.at(i).dst_to_origin = dst_to;
            if (dst_to < 5) { //原点の基地局にアクセスするデバイスを探す
                for (int j = 0; j < num_BS; j++) {
                    double dst = cal_dst(pos, BS_pos.at(j));
                    device.at(i).state = true;
                    if (dst < device.at(i).dst_to_origin) {
                        device.at(i).state = false;
                        break;
                    }
                    if (j == num_BS - 1) acl.push_back(i);
                }
            }
            //チャネル条件の逆数をかけて消えるので今は無視
            if (device.at(i).state) { //原点にアクセスするならallocateされた電力を足す．
                int level = L * device.at(i).dst_to_origin;;
                if (level == 0) level = 1;
                else if (level >= L) level = L;
                device.at(i).level = level;
                device.at(i).coef = theta * pow(theta + 1, L - level);
                LEV.at(level-1) += theta * pow(theta + 1, L - level);
            } else {
                //端末-基地局間のフェージング係数を設定
                double H = exp_dist(1.0);//gauss_rand(0, 1);
                double coef = H / pow(device.at(i).dst_to_origin, alpha);
                device.at(i).coef = coef;
                SI += coef;
            }
            
        }
        
        int s = (int)acl.size();
        for (int i = 0; i < s; i++) {
            int f = acl.at(i);
            double P = device.at(f).coef;
            int l = device.at(f).level;
            double si = 0;
            for (int j = 0; j < l; j++) {
                si += LEV.at(j);
            }
            double SINR = P / (SI + si - P + noise);
            if (SINR > theta) success++;
        }
        
        
        //進捗状況を表示
        double progress = i / end_time;
        if (progress > 0.8 && p_flag[3]) {cout << "...80%"; p_flag[3] = false;}
        else if (progress > 0.6 && p_flag[2]) {cout << "...60%"; p_flag[2] = false;}
        else if (progress > 0.4 && p_flag[1]) {cout << "...40%"; p_flag[1] = false;}
        else if (progress > 0.2 && p_flag[0]) {cout << "...20%"; p_flag[0] = false;}
    }
    cout << "...100%" << endl;
    
    double res = success / end_time;
    //double I = thp_theo(lambda_IoT, lambda_BS, alpha);
    cout << "Throughput : " <<  res << endl << endl;
    outputfile << res << " ";
}


//void decode_dst() {
//    //座標読み込み
//    string filename1 = "2.0Second_decode_pos.txt";
//    vector<pair<double, double>> BS_pos;
//    ifstream readingfile;
//    readingfile.open(filename1);
//    ifstream ifs("2.0Second_decode_pos.txt");
//        if (ifs.fail()) {
//           cerr << "Cannot open file\n";
//           exit(0);
//        }
//    ofstream outputfile3;
//    string filename3 = "2.0_4.5_secondm.txt";
//    outputfile3.open(filename3);
//    double x, y;
//    while (ifs >> x >> y) {
//        if (x <= 2.0 && y <= 2.0) {
//            outputfile3 << x << " " << y << endl;
//        }
//    }
//    ifs.close();
//
//}


int main() {
    cout << "Pos 0, Thp 1, SIC 2 :";
    double key; cin >> key; cout << endl;
    cout << "Put start lambda IoT : ";
    double k; cin >> k; cout << endl;
    
    string filename;
    if (key == 0) filename = to_string(k) + "PA_NOMA_pos.txt";
    else if (key == 1) filename = to_string(k) + "PA_NOMA_thp.txt";
    else {
        string filename2 = "Second_decode_pos.txt";
        string filename1 = "First_decode_pos.txt";
        outputfile2.open(filename2);
        outputfile1.open(filename1);
    }
    
    outputfile.open(filename);
    
//    for (double alpha = 2.5; alpha <= 4.5; alpha += 0.5) {
//        SIC_NOMA_thp(k, alpha, 0);
//    }
    
    for (double L = 3.0; L <= 5; L++) {
        double alpha = 4.5;
        if (key == 0) PA_NOMA_pos(k, alpha, 0, L);
        else PA_NOMA_thp(k, alpha, 0, L);
    }
    
    
    outputfile1.close();
    outputfile2.close();
    outputfile.close();
}
