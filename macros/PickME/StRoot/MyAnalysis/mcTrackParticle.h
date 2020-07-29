#ifndef MCTRACKPARTICLE_H
#define MCTRACKPARTICLE_H


#include "TMath.h"
#include "TObject.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "TVector2.h"

//class StMCTrackEvent;

class mcTrackParticle : public TObject {

 protected:
    // Track properties
    Double_t    Particle_pt;
    Double_t    Particle_eta;
    Double_t    Particle_phi;

    Double_t    dca_to_prim;
    Double_t    Particle_charge;

 public:

    mcTrackParticle();
    ~mcTrackParticle();

    void set_dca_to_prim(Double_t f)                     { dca_to_prim = f;            }
    void set_Particle_charge(Double_t f)                 { Particle_charge = f;        }
    //void set_TLV_Particle_prim(TLorentzVector tlv)      { TLV_Particle_prim = tlv;    }
    void set_Particle_pt(Double_t f)                     { Particle_pt = f;        }
    void set_Particle_eta(Double_t f)                    { Particle_eta = f;        }
    void set_Particle_phi(Double_t f)                    { Particle_phi = f;        }


    Double_t get_dca_to_prim()                      { return dca_to_prim;         }
    Double_t get_Particle_charge ()                 { return Particle_charge;         }
    //TLorentzVector get_TLV_Particle_prim()         { return TLV_Particle_prim;   }
    Double_t get_Particle_pt ()                      { return Particle_pt;         }
    Double_t get_Particle_eta ()                     { return Particle_eta;         }
    Double_t get_Particle_phi ()                     { return Particle_phi;         }

    ClassDef(mcTrackParticle,1)  // A simple track of a particle



};

#endif