#include <cstdlib>
#include <unistd.h>
#include <sstream>
#include <sys/types.h>
#include <sys/wait.h>
#include <regex>
#include "gateToFermat.h"

void string_replace_all(string &str, const string &from, const string &to) {
    size_t start_pos = 0;
    while((start_pos = str.find(from, start_pos)) != string::npos) {
        str.replace(start_pos, from.length(), to);
        start_pos += to.length();
    }
}

void string_trim(string &ostr) {
    const char* WhiteSpace = " \t\v\r\n";
    if(!ostr.empty()) {
        ostr.erase(0, ostr.find_first_not_of(WhiteSpace));
        ostr.erase(ostr.find_last_not_of(WhiteSpace)+1);
    }
}

class Fermat {
public:
    static int buffer_size;
    string Sentinel = "---EOF---";
    void Init(string fer_path="fer64");
    int vmax = 0;
    string Execute(string);
    void Exit();
    ~Fermat();
    
private:
    bool inited = false;
    bool exited = false;
    int P2C[2];
    int C2P[2];
    pid_t fpid = 0;
    pid_t pid = 0;
};

int Fermat::buffer_size = 1024*128;

#define ENTER endl<<endl<<endl

Fermat::~Fermat() { Exit(); }
    
void Fermat::Exit() {
    if(getpid()!=pid) return;
    if(inited) {
        ostringstream script;
        script << "&q;" << endl << "&x;" << ENTER;
        string istr = script.str();
        write(P2C[1], istr.c_str(), istr.length());
        int st;
        waitpid(fpid, &st, WUNTRACED);
        inited = false;
        exited = true;
    }
}

void Fermat::Init(string fer_path) {
    if(inited) return;
    inited = true;
    pid = getpid();
    
    if (pipe(P2C)==-1 || pipe(C2P)==-1) {
        cout << "pipe failed in Fermat::Init." << endl;
        abort();
    }
    
    fpid = fork();
    if (fpid == 0) { // child process
        if(setpgid(0,0)) {
            cout << "setpgid failed in Fermat::Init." << endl;
            abort();
        }
        close(P2C[1]);
        close(C2P[0]);
        dup2(C2P[1], 1);
        close(C2P[1]);
        dup2(P2C[0], 0);
        close(P2C[0]);
        execlp(fer_path.c_str(), fer_path.c_str(), NULL);
        exit(0);
    }
    
    // parent process
    close(P2C[0]);  // P2C[1] for write
    close(C2P[1]); // C2P[0] for read
    
    ostringstream script;
    script << "&(_d=90000);" << endl; // width of the display on the window
    script << "&d" << endl << "0;" << endl; // off floating point representation
    script << "&(_t=0);" << endl; // off a certain fast probabalistic algorithm
    script << "&(t=0);" << endl; // off timing
    script << "&(_s=0);" << endl;
    script << "&(_o=1000);" << endl; // http://home.bway.net/lewis/fer64mono.html
    script << "&(M=' ');" << endl; // prompt
    script << "!('" << Sentinel << "');" << ENTER;
    string istr = script.str();
    write(P2C[1], istr.c_str(), istr.length());
    
    string ostr;
    int n = 1024;
    char buffer[n+1]; // make sure the last one is '\0'
    int nio;
    while(true) {
        for(int i=0; i<n+1; i++) buffer[i] = '\0';
        nio = read(C2P[0], buffer, n);
        if(nio>0) ostr += buffer;
        else {
            cout << ostr << endl;
            abort();
        }
        auto cpos = ostr.find(Sentinel);
        if(cpos!=string::npos) {
            const char* WhiteSpace = " \t\v\r\n";
            auto lpos = ostr.find_last_not_of(WhiteSpace);
            if(ostr[lpos]!='0') read(C2P[0], buffer, n); // last 0, due to Sentinel
            ostr.erase(cpos);
            break;
        }
    }
    
    string_replace_all(ostr, "*** entry > 30 or < 5 means turn off mono multiply.", "");
    string_replace_all(ostr, "*** Fermat Warning. Early exit from mod_multivar_Chinese.", "");
    if(ostr.find("***")!=string::npos) {
        cout << "Fermat script: " << endl << istr << endl << endl;
        cout << ostr.c_str() << endl;
        abort();
    }
}

// out string still contains the last number
string Fermat::Execute(string expr) {
    if(exited) {
        cout << "Fermat has already exited." << endl;
        abort();
    }
    if(getpid() != pid) {
        cout << "Fermat: can not Execute on child process." << endl;
        abort();
    }
    ostringstream script;
    script << expr << endl;
    script << "!('" << Sentinel << "')" << ENTER;
    string istr = script.str();
    write(P2C[1], istr.c_str(), istr.length());
    
    string ostr;
    int n = buffer_size;
    char buffer[n+1]; // make sure the last one is '\0'
    int nio;
    while(true) {
        for(int i=0; i<n+1; i++) buffer[i] = '\0';
        nio = read(C2P[0], buffer, n);
        if(nio>0) ostr += buffer;
        else {
            cout << "input: " << expr << endl;
            cout << ostr << endl;
            abort();
        }
        auto cpos = ostr.find(Sentinel);
        if(cpos!=string::npos) {
            const char* WhiteSpace = " \t\v\r\n";
            auto lpos = ostr.find_last_not_of(WhiteSpace);
            if(ostr[lpos]!='0') read(C2P[0], buffer, n); // last 0, due to Sentinel
            ostr.erase(cpos);
            break;
        }
    }
    string_replace_all(ostr, "`", "");
    string_trim(ostr);

    if(ostr.find("***")!=string::npos) {
        cout << endl << expr << endl << endl;
        cout << ostr.c_str() << endl;
        abort();
    }
    return ostr;
}

static Fermat fermat;
static vector<string> v_vec;
static bool has_vv;
int iniCalc(const char *path, const string & vars) {
    fermat.Init(path);
    int n = std::count(vars.begin(), vars.end(), '\n');
    char cvs[vars.length()+1];
    strcpy(cvs, vars.c_str());
    char *pos = cvs;
    char *end = pos;
    vector<string> fvar_vec;
    has_vv = false;
    for(int i=0; i<n; i++) {
        while(*end != '\n') ++end;
        *end = '\0';
        string svar = string(pos);
        if(svar.find("fVAR")==0) fvar_vec.push_back(svar);
        if(std::isupper(svar[0])) {
            has_vv = true;
            v_vec.push_back(svar);
            svar = "fVAR"+svar;
        }
        fermat.Execute("&(J="+svar+");");
        pos = end;
        pos++;
    }
    if(has_vv && fvar_vec.size()>0) {
        cout << "Variables: ";
        for(auto vi : fvar_vec) cout << vi << " ";
        cout << endl;
        cout << "The variable name begins with fVAR, and it is not allowed!" << endl;
        abort();
    }
    return 0;
}

string calc(const string & in_buf) {
    string buf = in_buf;
    if(has_vv) {
        for(auto vi : v_vec) {
            std::regex vi_word("\\b" + vi + "\\b");
            buf = std::regex_replace(buf, vi_word, "fVAR"+vi);
        }
    }
    ostringstream oss;
    oss << "res:=" << buf << ";" << endl;
    oss << "&(U=1);" << endl; // ugly printing, the whitespace matters
    fermat.Execute(oss.str());
    oss.clear();
    oss.str("");
    oss << "!(^res):" << endl;
    auto ostr = fermat.Execute(oss.str());
    oss.clear();
    oss.str("");
    oss << "&(U=0);" << endl; // disable ugly printing
    oss << "@(res);" << endl;
    oss << "&_G;" << endl;
    fermat.Execute(oss.str());
    if(has_vv) {
        for(auto vi : v_vec) {
            std::regex fvi_word("\\bfVAR" + vi + "\\b");
            ostr = std::regex_replace(ostr, fvi_word, vi);
        }
    }
    return ostr;
}

void closeCalc() {
    fermat.Exit();
}
