#include <iostream>
#include <fstream>
#include <regex>
#include <vector>
#include <map>
#include <algorithm>

// outputs both

class barrierAOE {
public:
    int getStart() {
        return m_start;
    }
    
    barrierAOE(int start, int zero, int stop, char orientation) {
        m_start = start;
        m_zero = zero;
        m_stop = stop;
        m_orientation = orientation;
    }
    
    int isInside(int pos) {
        if(pos >= m_start && pos <= m_stop) {
            return 0;
        } else if (pos < m_start) {
            return 1;
        } else if (pos > m_stop) {
            return 2;
        } else {
            std::cout<<"Start :"<<m_start<<" & stop : "<<m_stop<<std::endl;
            throw std::domain_error("m_stop < m_start");
        }
    }
    
    int distToBorder(int pos) {
        if(m_orientation == 'R') {
            pos = -1 * pos;
        }
        return pos - m_zero;
    }
    
private:
    int m_start;
    int m_zero;
    int m_stop;
    int m_orientation;
};

bool sortAOE(barrierAOE i, barrierAOE j) {
    return i.getStart() < j.getStart();
}

int main(int argc, char* argv[]) {
    std::ifstream AOEFile(argv[1]);
    std::ifstream mutedFile(argv[2]);
    std::ofstream countFile(argv[3]);
    
    char tchar = 'X';
    std::string chrom = "";
    std::string startS = "";
    std::string zeroS = "";
    std::string stopS = "";
    char orientation = '+';
    int start = 0;
    
    std::map <int, std::map <char, int>> basesByPos;
    
    std::map <std::string, std::vector <barrierAOE>> barriersAOE;
    int col = 0;
    
    while(!AOEFile.eof()) {
        AOEFile.get(tchar);
        // file is chrom : start : zero : stop
        if(tchar == '\t') {
            col ++;
        } else if (tchar != '\n') {
            switch(col) {
                case 0:
                    chrom += tchar;
                    break;
                case 1:
                    startS += tchar;
                    break;
                case 2:
                    zeroS += tchar;
                    break;
                case 3:
                    stopS += tchar;
                    break;
                case 4:
                    orientation = tchar;
                    break;
                default:
                    std::cout<<col<<std::endl;
                    throw std::domain_error("unexpected value for col var");
            } 
        } else if(!AOEFile.eof()) { // necessarly a '\n'
            // std::cout << startS << " " << zeroS << " " << stopS << " " << orientation << std::endl;
            barrierAOE temp(stoi(startS), stoi(zeroS), stoi(stopS), orientation);
            if(barriersAOE.find(chrom) == barriersAOE.end())  {
                barriersAOE[chrom] = std::vector <barrierAOE>(1, temp);
            } else {
                barriersAOE[chrom].push_back(temp);
            }
            chrom = "";
            startS = "";
            zeroS = "";
            stopS = "";
            col = 0;
        }
    }
    
    for(const auto &pair : barriersAOE) {
        std::string key = pair.first;
        std::sort(barriersAOE[key].begin(), barriersAOE[key].end(), sortAOE);
        // std::cout<< barriersAOE[key].size() << std::endl;
    }
    
    tchar = 'X';
    bool isStart = false;
    int line = 0;
    
    while(!mutedFile.eof()) {
        mutedFile.get(tchar);
        if(tchar == '>') {
            isStart = false;
            chrom = "";
            startS = "";
            while(tchar != '\n') {
                mutedFile.get(tchar);
                if(tchar == ':') {
                    isStart = true;
                } else if(isStart) {
                    startS += tchar;
                } else {
                    chrom += tchar;
                }
            }
            // std::cout<<chrom<< " " << startS<<std::endl;
            start = stoi(startS);
            
        } else if(tchar != '\n') {
            start ++;
            if(tchar != '*' && tchar != '-') {
                // need to compare to array of barriers AOE
                //std::vector <barrierAOE> array = barriersAOE[chrom];
                int startI = 0;
                int stopI = barriersAOE[chrom].size();
                bool found = false;
                
                while(stopI - startI > 2 && !found) {
                    int middle = (startI+stopI) /2;
                    int inside = barriersAOE[chrom][middle].isInside(start);
                    switch(inside) {
                        case 0:
                            basesByPos[barriersAOE[chrom][middle].distToBorder(start)][toupper(tchar)] ++;
                            found = true;
                            break;
                        case 1: // smaller
                            stopI = middle;
                            break;
                        case 2: // bigger
                            startI = middle;
                            break;
                        default:
                            std::cout<<"WTF"<<std::endl;
                            break;
                    }
                }
            }
        } else {
            line++;
            if(line % 1000 == 0) {
                std::cout<<line<<"                           \r";
            }
        }
    }
    std::cout<<std::endl;
    countFile << "pos\tA\tB\tC\tD\tE\tF\tG\tH\tI\tJ\tK\tL\tM\tN\tO\tP\tQ\tR\tS\tT\tV" << std::endl;
//     for(const auto &pair : basesByPos) {
//         for(const auto &pair2 : pair.second) {
//             std::cout<<pair.first<<"  bases : "<<pair2.first<<std::endl;
//         }
//     }
    
    for(const auto &pair : basesByPos) {
        std::map <char, int> pos = pair.second;
        std::string bases = "ABCDEFGHIJKLMNOPQRSTV";
        countFile << pair.first;
        for(int i(0); i< bases.size(); i++) {
            countFile << '\t' << pos[bases[i]];
        }
        countFile << '\n';
    }
    
    return 0;
}
