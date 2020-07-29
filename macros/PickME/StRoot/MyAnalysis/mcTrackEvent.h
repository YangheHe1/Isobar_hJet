#ifndef MCTRACKEVENT_H
#define MCTRACKEVENT_H

#include "TVector3.h"
#include "TClonesArray.h"
#include "TObject.h"
#include <vector>
#include "mcTrackParticle.h"

class mcTrackParticle;

class mcTrackEvent : public TObject {
 protected:

    Double_t x;
    Double_t y;
    Double_t z;
    Int_t   id;
    Double_t mult;
    Double_t numtrack;
    Double_t psi;
    
    Double_t ZDCx;
    Double_t BBCx;
    Double_t vzVpd;

    UShort_t      fNumParticle;

    TClonesArray* fParticle;

 public:

    mcTrackEvent();
    ~mcTrackEvent();

        void       setx(Double_t r)                    { x = r;                         }
        Double_t    getx()                        { return x;                      }

        void       sety(Double_t r)                    { y = r;                         }
        Double_t    gety()                        { return y;                      }

        void       setz(Double_t r)                    { z = r;                         }
        Double_t    getz()                        { return z;                      }

        void       setid(Int_t  r)                    { id = r;                        }
        Int_t      getid()                       { return id;                     }

        void       setmult(Double_t r)                 { mult = r;                      }
        Double_t    getmult()                     { return mult;                   }

        void       setnumtrack(Double_t r)                 { numtrack = r;                      }
        Double_t    getnumtrack()                     { return numtrack;                   }

        void       setpsi(Double_t r)                 { psi = r;                      }
        Double_t    getpsi()                     { return psi;                   }

        void       setZDCx(Double_t r)                 { ZDCx = r;                      }
        Double_t    getZDCx()                     { return ZDCx;                   }

        void       setBBCx(Double_t r)                 { BBCx = r;                      }
        Double_t    getBBCx()                     { return BBCx;                   }

        void       setvzVpd(Double_t r)                { vzVpd = r;                     }
        Double_t    getvzVpd()                    { return vzVpd;                  }

        mcTrackParticle* createParticle();
    
        void clearParticleList()
        {
            fNumParticle   = 0;
            fParticle      ->Clear();
        }

         UShort_t getNumParticle() 
        {
            return fNumParticle;
        }

        mcTrackParticle* getParticle(UShort_t i) 
        {
            return i < fNumParticle ? (mcTrackParticle*)((*fParticle)[i]) : NULL;
        }
        //--------------------------------------


        ClassDef(mcTrackEvent,1)  // A simple event compiled of tracks




};

#endif