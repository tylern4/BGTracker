////////////////////////////////////////////////////////////////////
//
// Dean Andrew Hidas <dhidas@bnl.gov>
//
// Created on: Tue Aug 23 11:37:14 EDT 2022
//
////////////////////////////////////////////////////////////////////


#include <iostream>
#include <vector>

#include "Drif.h"
#include "Quad.h"
#include "Bend.h"
#include "Sext.h"
#include "BGCuda.h"

std::vector<BGE*> BuildLattice ();



int TrackingTest (int const NTurns, int const NParticles)
{
  // oh.. forgive me for this mess

  int const NThreadsPerBlock = 192;
  int const NBlocks = NParticles % NThreadsPerBlock == 0 ? ((int) NParticles / NThreadsPerBlock) : ((int) NParticles / NThreadsPerBlock) + 1;
  std::cout << "NBlocks: " << NBlocks << std::endl;
  std::cout << NParticles*6*sizeof(double) / 1024/1024 << " MB particle data"  << std::endl;

  // Set a random seed
  std::srand(1234987);

  // Get particle data
  double* particles;
  cudaMallocManaged(&particles, NParticles*6*sizeof(double));

  // Can randomize this later..
  for (int i = 0; i != NParticles; ++i) {
      particles[i*6 + 0] = 1e-8;
      particles[i*6 + 1] = 0;
      particles[i*6 + 2] = 0;
      particles[i*6 + 3] = 0;
      particles[i*6 + 4] = 0;
      particles[i*6 + 5] = 0;
      //for (int j = 0; j != 6; ++j) {
      //  particles[i*6 + j] = (double)std::rand()/(float) RAND_MAX;
      //}
  }


  std::vector<BGE*> Lattice = BuildLattice();
  std::cout << "Number of elements not including NKick for each: " << Lattice.size() << std::endl;

  cudaDeviceSynchronize();

  for (int iturn = 0; iturn != NTurns; ++iturn) {
    for (std::vector<BGE*>::iterator it = Lattice.begin(); it != Lattice.end(); ++it) {
      switch ((*it)->Type()) {
        case BGE::ElementType::Drif:
          {
            Drif* d = (Drif*) *it;
            dsympass4<<<NBlocks, NThreadsPerBlock>>>(particles, d->fTM, d->fL, NParticles);
            break;
          }
        case BGE::ElementType::Quad:
          {
            Quad* q = (Quad*) *it;
            qsympass4<<<NBlocks, NThreadsPerBlock>>>(particles, q->fMa, q->fMb, q->fK1Lg, q->fK1Ld, q->fdL, q->fNKick, NParticles);
            break;
          }
        case BGE::ElementType::Bend:
          {
            Bend* b = (Bend*) *it;
            bsympass4<<<NBlocks, NThreadsPerBlock>>>(particles, b->fMa, b->fMb, b->fM1, b->fM2, b->fK1Ld, b->fK1Lg, b->fK2Ld, b->fK2Lg, b->fLg, b->fLd,  b->fdL, b->fR, b->fNKick, NParticles);
            break;
          }
        case BGE::ElementType::Sext:
          {
            Sext* s = (Sext*) *it;
            ssympass4<<<NBlocks, NThreadsPerBlock>>> (particles, s->fMa, s->fMb, s->fK2Ld, s->fK2Lg, s->fdL, s->fNKick, NParticles);
            break;
          }
      }
    }
  }

  cudaDeviceSynchronize();


  for (int i = 0; i != 1; ++i) {
    for (int j = 0; j != 6; ++j) {
      printf("%.6e ", particles[i*6 + j]);
    }
    std::cout << std::endl;
  }

  cudaFree(particles);
  return 0;
}


int main (int argc, char* argv[])
{
  if (argc != 3) {
    std::cerr << "Usage: " << argv[0] << " [NTurns] [NParticles]" << std::endl;
    return 1;
  }

  int const NTurns = atoi(argv[1]);
  int const NParticles = atoi(argv[2]);
  TrackingTest(NTurns, NParticles);


  return 0;
}





std::vector<BGE*> BuildLattice ()
{
  std::vector<BGE*> L;
  Bend *b1g3c01a = new Bend("b1g3c01a", 2.62, 0.10471975512, 0.05236, 0.05236, 0.0, 0.0);
  Bend *b1g3c30a = new Bend("b1g3c30a", 2.62, 0.10471975512, 0.05236, 0.05236, 0.0, 0.0);
  Bend *b1g5c01b = new Bend("b1g5c01b", 2.62, 0.10471975512, 0.05236, 0.05236, 0.0, 0.0);
  Bend *b1g5c30b = new Bend("b1g5c30b", 2.62, 0.10471975512, 0.05236, 0.05236, 0.0, 0.0);
  Drif *D0001 = new Drif("D0001", 4.65);
  Drif *D0002 = new Drif("D0002", 0.166);
  Drif *D0003 = new Drif("D0003", 0.802);
  Drif *D0004 = new Drif("D0004", 0.184);
  Drif *D0005 = new Drif("D0005", 0.186);
  Drif *D0006 = new Drif("D0006", 0.166);
  Drif *D0007 = new Drif("D0007", 0.58);
  Drif *D0008 = new Drif("D0008", 0.9765);
  Drif *D0009 = new Drif("D0009", 0.2015);
  Drif *D0010 = new Drif("D0010", 0.654);
  Drif *D0011 = new Drif("D0011", 0.184);
  Drif *D0012 = new Drif("D0012", 0.184);
  Drif *D0013 = new Drif("D0013", 0.504);
  Drif *D0014 = new Drif("D0014", 0.3515);
  Drif *D0015 = new Drif("D0015", 0.9765);
  Drif *D0016 = new Drif("D0016", 0.591);
  Drif *D0017 = new Drif("D0017", 0.166);
  Drif *D0018 = new Drif("D0018", 0.5596);
  Drif *D0019 = new Drif("D0019", 0.244);
  Drif *D0020 = new Drif("D0020", 0.466);
  Drif *D0021 = new Drif("D0021", 0.166);
  Drif *D0022 = new Drif("D0022", 6.6);
  Drif *D0023 = new Drif("D0023", 0.166);
  Drif *D0024 = new Drif("D0024", 0.466);
  Drif *D0025 = new Drif("D0025", 0.244);
  Drif *D0026 = new Drif("D0026", 0.5596);
  Drif *D0027 = new Drif("D0027", 0.166);
  Drif *D0028 = new Drif("D0028", 0.591);
  Drif *D0029 = new Drif("D0029", 0.9765);
  Drif *D0030 = new Drif("D0030", 0.2015);
  Drif *D0031 = new Drif("D0031", 0.654);
  Drif *D0032 = new Drif("D0032", 0.184);
  Drif *D0033 = new Drif("D0033", 0.184);
  Drif *D0034 = new Drif("D0034", 0.504);
  Drif *D0035 = new Drif("D0035", 0.3515);
  Drif *D0036 = new Drif("D0036", 0.9765);
  Drif *D0037 = new Drif("D0037", 0.58);
  Drif *D0038 = new Drif("D0038", 0.166);
  Drif *D0039 = new Drif("D0039", 0.186);
  Drif *D0040 = new Drif("D0040", 0.184);
  Drif *D0041 = new Drif("D0041", 0.802);
  Drif *D0042 = new Drif("D0042", 0.166);
  Drif *D0043 = new Drif("D0043", 4.65);
  Quad *qh1g2c30a = new Quad("qh1g2c30a", 0.268, -0.641957314648);
  Quad *qh1g6c01b = new Quad("qh1g6c01b", 0.268, -0.641957314648);
  Quad *qh2g2c30a = new Quad("qh2g2c30a", 0.46, 1.43673057073);
  Quad *qh2g6c01b = new Quad("qh2g6c01b", 0.46, 1.43673057073);
  Quad *qh3g2c30a = new Quad("qh3g2c30a", 0.268, -1.75355042529);
  Quad *qh3g6c01b = new Quad("qh3g6c01b", 0.268, -1.75355042529);
  Quad *ql1g2c01a = new Quad("ql1g2c01a", 0.268, -1.61785473561);
  Quad *ql1g6c30b = new Quad("ql1g6c30b", 0.268, -1.61785473561);
  Quad *ql2g2c01a = new Quad("ql2g2c01a", 0.46, 1.76477357129);
  Quad *ql2g6c30b = new Quad("ql2g6c30b", 0.46, 1.76477357129);
  Quad *ql3g2c01a = new Quad("ql3g2c01a", 0.268, -1.51868267756);
  Quad *ql3g6c30b = new Quad("ql3g6c30b", 0.268, -1.51868267756);
  Quad *qm1g4c01a = new Quad("qm1g4c01a", 0.247, -0.812234822773);
  Quad *qm1g4c01b = new Quad("qm1g4c01b", 0.247, -0.812234822773);
  Quad *qm1g4c30a = new Quad("qm1g4c30a", 0.247, -0.812234822773);
  Quad *qm1g4c30b = new Quad("qm1g4c30b", 0.247, -0.812234822773);
  Quad *qm2g4c01a = new Quad("qm2g4c01a", 0.282, 1.22615465959);
  Quad *qm2g4c01b = new Quad("qm2g4c01b", 0.282, 1.22615465959);
  Quad *qm2g4c30a = new Quad("qm2g4c30a", 0.282, 1.22615465959);
  Quad *qm2g4c30b = new Quad("qm2g4c30b", 0.282, 1.22615465959);
  Sext *sh1  = new Sext("sh1", 0.2, 19.8329120997);
  Sext *sh3  = new Sext("sh3", 0.2, -5.85510841147);
  Sext *sh4  = new Sext("sh4", 0.2, -15.8209007067);
  Sext *sl1  = new Sext("sl1", 0.2, -13.2716060547);
  Sext *sl2  = new Sext("sl2", 0.2, 35.6779214531);
  Sext *sl3  = new Sext("sl3", 0.2, -29.4608606061);
  Sext *sm1a = new Sext("sm1a", 0.2, -23.6806342393);
  Sext *sm1b = new Sext("sm1b", 0.2, -25.9460354618);
  Sext *sm2  = new Sext("sm2", 0.25, 28.6431546915);
  //rfc = rfca("rfc",L=0,voltage=3e6,freq=0.49968e9)


  for (int i = 0; i != 15; ++i) {
    L.push_back( (BGE*) D0001);
    L.push_back( (BGE*) sh1);
    L.push_back( (BGE*) D0002);
    L.push_back( (BGE*) qh1g2c30a);
    L.push_back( (BGE*) D0003);
    L.push_back( (BGE*) qh2g2c30a);
    L.push_back( (BGE*) D0004);
    L.push_back( (BGE*) sh3);
    L.push_back( (BGE*) D0005);
    L.push_back( (BGE*) qh3g2c30a);
    L.push_back( (BGE*) D0006);
    L.push_back( (BGE*) sh4);
    L.push_back( (BGE*) D0007);
    L.push_back( (BGE*) b1g3c30a);
    L.push_back( (BGE*) D0008);
    L.push_back( (BGE*) qm1g4c30a);
    L.push_back( (BGE*) D0009);
    L.push_back( (BGE*) sm1a);
    L.push_back( (BGE*) D0010);
    L.push_back( (BGE*) qm2g4c30a);
    L.push_back( (BGE*) D0011);
    L.push_back( (BGE*) sm2);
    L.push_back( (BGE*) D0012);
    L.push_back( (BGE*) qm2g4c30b);
    L.push_back( (BGE*) D0013);
    L.push_back( (BGE*) sm1b);
    L.push_back( (BGE*) D0014);
    L.push_back( (BGE*) qm1g4c30b);
    L.push_back( (BGE*) D0015);
    L.push_back( (BGE*) b1g5c30b);
    L.push_back( (BGE*) D0016);
    L.push_back( (BGE*) ql3g6c30b);
    L.push_back( (BGE*) D0017);
    L.push_back( (BGE*) sl3);
    L.push_back( (BGE*) D0018);
    L.push_back( (BGE*) ql2g6c30b);
    L.push_back( (BGE*) D0019);
    L.push_back( (BGE*) sl2);
    L.push_back( (BGE*) D0020);
    L.push_back( (BGE*) ql1g6c30b);
    L.push_back( (BGE*) D0021);
    L.push_back( (BGE*) sl1);
    L.push_back( (BGE*) D0022);
    L.push_back( (BGE*) sl1);
    L.push_back( (BGE*) D0023);
    L.push_back( (BGE*) ql1g2c01a);
    L.push_back( (BGE*) D0024);
    L.push_back( (BGE*) sl2);
    L.push_back( (BGE*) D0025);
    L.push_back( (BGE*) ql2g2c01a);
    L.push_back( (BGE*) D0026);
    L.push_back( (BGE*) sl3);
    L.push_back( (BGE*) D0027);
    L.push_back( (BGE*) ql3g2c01a);
    L.push_back( (BGE*) D0028);
    L.push_back( (BGE*) b1g3c01a);
    L.push_back( (BGE*) D0029);
    L.push_back( (BGE*) qm1g4c01a);
    L.push_back( (BGE*) D0030);
    L.push_back( (BGE*) sm1a);
    L.push_back( (BGE*) D0031);
    L.push_back( (BGE*) qm2g4c01a);
    L.push_back( (BGE*) D0032);
    L.push_back( (BGE*) sm2);
    L.push_back( (BGE*) D0033);
    L.push_back( (BGE*) qm2g4c01b);
    L.push_back( (BGE*) D0034);
    L.push_back( (BGE*) sm1b);
    L.push_back( (BGE*) D0035);
    L.push_back( (BGE*) qm1g4c01b);
    L.push_back( (BGE*) D0036);
    L.push_back( (BGE*) b1g5c01b);
    L.push_back( (BGE*) D0037);
    L.push_back( (BGE*) sh4);
    L.push_back( (BGE*) D0038);
    L.push_back( (BGE*) qh3g6c01b);
    L.push_back( (BGE*) D0039);
    L.push_back( (BGE*) sh3);
    L.push_back( (BGE*) D0040);
    L.push_back( (BGE*) qh2g6c01b);
    L.push_back( (BGE*) D0041);
    L.push_back( (BGE*) qh1g6c01b);
    L.push_back( (BGE*) D0042);
    L.push_back( (BGE*) sh1);
    L.push_back( (BGE*) D0043);
  }


  return L;
}
