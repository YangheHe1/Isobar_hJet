#include <vector>

void paixu(){

	vector<int> shuzu;
	std::ifstream nfile ;
        nfile.open("runnumber.list");

        std::string nline;
        while (std::getline(nfile, nline))
        {
                std::istringstream iss(nline);
                int nf;
                vector<int> p__;
                while (iss >> nf)
                {
                        if(nf>0)        p__.push_back(nf);
                }


                if(p__.size()<=0) continue;
                std::cout <<p__.at(0)<<"\n";
                shuzu.push_back(p__.at(0));
                vector<int>().swap(p__);
        }
        nfile.close();

for(int m=0;m<shuzu.size()-1;m++){
for(int i=0;i<shuzu.size()-m-1;i++){
int zhongjian;
if(shuzu.at(i)>shuzu.at(i+1)){zhongjian=shuzu.at(i+1);shuzu.at(i+1)=shuzu.at(i);shuzu.at(i)=zhongjian;}
}
}

	ofstream fout;
        fout.open("runnumber1.list");
	for(int l=0;l<shuzu.size();l++) fout<<shuzu.at(l)<<endl;
	fout.close();

}
