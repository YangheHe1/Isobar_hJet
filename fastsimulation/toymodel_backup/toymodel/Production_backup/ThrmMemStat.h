/***************************************************************************
 *
 * $Id: ThrmMemStat.h,v 1.1 2009/01/26 14:31:48 fisyak Exp $
 *
 * Author: Victor Perev, Jul 2000
 **************************************************************************/
#ifndef ThrmMemStat_h
#define ThrmMemStat_h
#include "TNamed.h"

class TList;

class ThrmMemStat :public TNamed {
public:
    ThrmMemStat(const char *name=0);
   ~ThrmMemStat();
   void Start();
   void Stop();
   virtual void   Print(const char *tit="") const;

   //static methods

   static  Double_t Used();                     
   static  Double_t Free();                     
   static  Double_t ProgSize();                 
   static  void     PrintMem(const char *tit);  
   static  void     PM();                       
   static  void     Summary();                  
 private:
   Double_t fLast;
   Double_t fMin;
   Double_t fAver;
   Double_t fMax;
   Double_t fRms;
   Int_t    fTally;

   static Double_t fgUsed;
   static TList    *fgList;

   ClassDef(ThrmMemStat,0)
};

#endif

