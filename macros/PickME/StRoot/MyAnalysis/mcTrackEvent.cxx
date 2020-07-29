#include "mcTrackEvent.h"
#include "stdio.h"

#include "mcTrackParticle.h"

ClassImp(mcTrackEvent)

mcTrackEvent::mcTrackEvent() : 
        x(-1),y(-1),z(-1),id(-1),mult(0),numtrack(0),psi(0),
        ZDCx(-1),BBCx(-1),vzVpd(-1),fNumParticle(0)
{
    fParticle      = new TClonesArray( "mcTrackParticle", 10 );
}

mcTrackEvent::~mcTrackEvent(){
    delete fParticle;
    fParticle = NULL;
}

mcTrackParticle* mcTrackEvent::createParticle()
{
            if (fNumParticle == fParticle->GetSize())
                fParticle->Expand( fNumParticle + 5 );
            if (fNumParticle >= 5000)
            {
                Fatal( "mcTrackParticle::createParticle()", "ERROR: Too many Particles (>5000)!" );
                exit( 2 );
            }

            new((*fParticle)[fNumParticle++]) mcTrackParticle;
            return (mcTrackParticle*)((*fParticle)[fNumParticle - 1]);
}
