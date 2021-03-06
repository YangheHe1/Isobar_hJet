/***************************************************************************
 *
 * $Id: ThrmMemStat.cxx,v 1.2 2009/08/28 16:38:55 fine Exp $
 *
 ***************************************************************************
 *
 * Description:
 *
 ***************************************************************************
 **************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <unistd.h>
#include "ThrmMemStat.h"
#include "TList.h"
#include "TError.h"
#include <cassert>

Double_t  ThrmMemStat::fgUsed=0;
TList    *ThrmMemStat::fgList=0;
ClassImp(ThrmMemStat)

#define LOWEST_VAL 0.0000001                                           
//______________________________________________________________________________
ThrmMemStat::ThrmMemStat(const char *name):TNamed(name,"")
{
  int n = (char*)&fTally - (char*)&fLast + sizeof(fTally);
  memset(&fLast,0,n);
  fMin =  1.e+33; 
  fMax = -1.e+33; 
  if (!fgList) fgList=new TList;
  fgList->Add(this);
}
//______________________________________________________________________________
ThrmMemStat::~ThrmMemStat()
{
  fgList->Remove(this);
  if (!fgList->First()) {delete fgList; fgList=0;}
}

//______________________________________________________________________________
void ThrmMemStat::Start()
{
  fLast = Used();
}

//______________________________________________________________________________
void ThrmMemStat::Stop()
{
  fTally++;
  Double_t dif = Used() - fLast;

  //printf("DEBUG >> time distance between two stops Used=%f Last=%f\n",Used(),fLast);
  if ( fabs(dif) < LOWEST_VAL )  dif  = 0.0;
  if ( dif < fMin )        fMin = dif;
  if ( dif > fMax )        fMax = dif;

  fAver += dif;
  fRms  += (dif*dif);

}
//______________________________________________________________________________
void ThrmMemStat::Print(const char *) const
{
  if (!fTally) return;
  Double_t aver = fAver/fTally;
  Double_t rms  = ::sqrt(fabs(fRms/fTally - aver*aver));

  //printf("DEBUG :: %.10f %d %.10f %.10f\n",fAver,fTally,fRms,aver);
  if ( fabs(aver) < LOWEST_VAL ) aver = 0.0;
  if ( rms        < LOWEST_VAL ) rms  = 0.0;

  printf("%40.40s(%d)%12.6f%12.6f%12.6f%12.6f\n",
         GetName(),fTally,fMin,aver,fMax,rms);
}
//______________________________________________________________________________
void ThrmMemStat::Summary()
{
#define NUMTICKS (40+4*12+5)

  Double_t dmin=1.e+33,daver=0,dmax=-1.e+33,drms=0,dtally=0,dmp;
  int i;

  if(!fgList) return;
  fgList->Sort();
  printf("%40.40s%12s%12s%12s%12s\n",
         "ThrmMemStat::Summary(calls)","Min ","Aver ","Max ","RMS ");

  for( i=0 ; i < NUMTICKS ; i++) printf("=");
  printf("\n");
   
  TListIter next(fgList); 
  ThrmMemStat  *m;
  while((m = (ThrmMemStat*)next())){
    if(!m->fTally)      continue;
    m->Print();
    dtally++;
    if (m->fMin < dmin) dmin=m->fMin;
    if (m->fMax > dmax) dmax=m->fMax;
    dmp = m->fAver/m->fTally;
    daver += dmp; 
    drms  += fabs(m->fRms/m->fTally-dmp*dmp);
    
  }
  if(!dtally) return;

  for( i=0 ; i < NUMTICKS ; i++) printf("-");
  printf("\n");

  //VP daver /= dtally;
  drms   = ::sqrt(fabs(drms));
  printf("%40.40s(%d)%12.6f%12.6f%12.6f%12.6f\n",
          "Total",(int)dtally,dmin,daver,dmax,drms);

  for( i=0 ; i < NUMTICKS ; i++) printf("=");
  printf("\n");

}

//______________________________________________________________________________
Double_t ThrmMemStat::Used()
{
  struct mallinfo info;
  info = mallinfo();
  return double(info.uordblks + info.usmblks)/1000000;
}
//______________________________________________________________________________
Double_t ThrmMemStat::Free()
{
  struct mallinfo info;
  info = mallinfo();
  return double(info.fordblks + info.fsmblks)/1000000;
}

//______________________________________________________________________________
Double_t ThrmMemStat::ProgSize()
{
  Double_t    res=0;
  int pid = ::getpid();
  char line[100];
  sprintf(line,"/proc/%d/status",pid);

FILE *proc = fopen(line,"r");
  if (proc) {//status file found
    while (fgets(line,100,proc)) {
      if (strncmp("VmSize:",line,7)==0) {
        fclose(proc); 
        char *aft=0;
        res = (strtod(line+7,&aft));  
        while ((++aft)[0]==' '){}
        int b = 0;
        if (strncmp("kB",aft,2)==0) b = 1024;
        if (strncmp("mB",aft,2)==0) b = 1024*1024;
        if (strncmp("gB",aft,2)==0) b = 1024*1024*1024;
        res = (res*b)/(1024*1024);
        return res;
      }
    }
    fclose(proc);
  }
//    status file not found. Use ugly way via "ps"
  static char *ps = 0;
  if (!ps) {
    ps = (char*)malloc(25);
    sprintf(ps,"/bin/ps -l -p %d",pid);
  }
  FILE *pipe = ::popen(ps,"r");
  if (!pipe) return 0.;
  
  char   psBuf[130];
  psBuf[0] = ' ';
  while( !feof( pipe ) ) {
    psBuf[1]=0;
    if(!fgets( psBuf+1, 128, pipe)) continue;
//    printf("%s\n",psBuf);
    int ifild=0;char *c;

    for (c=psBuf; c[0]; c++) {
      if (c[0]==' ' && c[1]!=' ') ifild++;
      if (ifild == 10) break;
    }
    res = (Double_t) atoi(c+1);
    if (res) break;
  }
  ::pclose(pipe);
  res *=::getpagesize()/(1024.*1024.);

  return res;
}

//______________________________________________________________________________
void ThrmMemStat::PrintMem(const char *tit)
{
  Double_t used = Used();
  Double_t free = Free();
  Double_t exec = ProgSize();

  if (tit) printf("\nThrmMemStat::%s",tit);
  printf("\t total =%10.6f heap =%10.6f and %10.6f(%+10.6f)\n",exec,used,free,used-fgUsed);
  fgUsed = used;
}
//______________________________________________________________________________
void ThrmMemStat::PM()
{
  Double_t used = Used();
  printf("\nThrmMemStat: ");
  printf("ThrmMemStat::heap =%10.6f(%+10.6f)\n",used,used-fgUsed);
  fgUsed = used;
}
//______________________________________________________________________________
//void ThrmMemStat::Streamer(TBuffer&)
//{assert(0);}
//______________________________________________________________________________
// void ThrmMemStat::ShowMembers(TMemberInspector& insp, char* parent){}
