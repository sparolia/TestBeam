
// Daniel Pitzl, DESY, Feb 2016, Jan 2018, May 2019
// telescope analysis with eudaq and openMP
// triplet alignment

// make tele
// tele -l 99999 33095
// needs runs.dat

#include "eudaq/FileReader.hh"
#include "eudaq/PluginManager.hh"
#include "../main/lib/plugins/BDAQ53ConverterPlugin.cc"

#include <TFile.h>
#include <TH1I.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TF1.h>

#include <sstream> // stringstream
#include <fstream> // filestream
#include <list>
#include <map>
#include <set>
#include <cmath>
#include <time.h> // clock_gettime

using namespace std;
using namespace eudaq;

struct pixel {
  int col;
  int row;
  int nn;
};

struct cluster {
  vector <pixel> vpix;
  //int size; int ncol, nrow;
  unsigned scr; // compressed
  float col, row;
  float mindxy;
};

struct triplet {
  double xm; // mid
  double ym;
  double zm;
  double sx; // slope
  double sy;
  double rxy; // residual
  double kx; // kink
  double ky;
  int ghost;
  unsigned i1; // cluster index
  unsigned i2;
  unsigned i3;
  vector <double> vx; // x cluster
  vector <double> vy;
};

struct vertex {
  double x;
  double y;
  double z;
  unsigned i;
  unsigned j;
};

bool ldbg = 0; // global

//------------------------------------------------------------------------------
vector<cluster> getClus( vector <pixel> pb, int fCluCut = 1 ) // 1 = no gap
{
  // returns clusters with pixel coordinates
  // next-neighbour topological clustering (allows fCluCut-1 empty pixels)

  vector <cluster> vc;
  if( pb.size() == 0 ) return vc;

  int * gone = new int[pb.size()];
  for( unsigned i = 0; i < pb.size(); ++i )
    gone[i] = 0;

  unsigned seed = 0;

  while( seed < pb.size() ) {

    // start a new cluster:

    cluster c;
    c.vpix.push_back( pb[seed] );
    gone[seed] = 1;

    // let it grow as much as possible:

    int growing;
    do{
      growing = 0;
      for( unsigned i = 0; i < pb.size(); ++i ) {
        if( !gone[i] ){ // unused pixel
          for( unsigned int p = 0; p < c.vpix.size(); ++p ) { // vpix in cluster so far
            int dr = c.vpix.at(p).row - pb[i].row;
            int dc = c.vpix.at(p).col - pb[i].col;
	    //if( (   dr>=-fCluCut) && (dr<=fCluCut) 
	    //&& (dc>=-fCluCut) && (dc<=fCluCut) ) { // allow diagonal
            if( ( abs(dr) <= fCluCut && dc == 0 ) ||
		( abs(dc) <= fCluCut && dr == 0 ) ) { // only facing neighbours, same resolution
              c.vpix.push_back(pb[i]);
	      gone[i] = 1;
              growing = 1;
              break; // important!
            }
          } // loop over vpix
        } // not gone
      } // loop over all pix
    }
    while( growing );

    // count pixel neighbours:

    for( vector <pixel>::iterator p = c.vpix.begin(); p != c.vpix.end(); ++p ) {
      vector <pixel>::iterator q = p;
      ++q;
      for( ; q != c.vpix.end(); ++q )
	if( abs( p->col - q->col ) <= 1 &&abs( p->row - q->row ) <= 1 ) {
	  ++p->nn;
	  ++q->nn;
	}
    }

    c.col = 0;
    c.row = 0;
    int sumnn = 0;
    int minx = 9999;
    int maxx = 0;
    int miny = 9999;
    int maxy = 0;

    for( vector <pixel>::iterator p = c.vpix.begin(); p != c.vpix.end(); ++p ) {

      int nn = max( 1, p->nn ); // neighbours
      sumnn += nn;
      c.col += p->col * nn;
      c.row += p->row * nn;
      if( p->col > maxx ) maxx = p->col;
      if( p->col < minx ) minx = p->col;
      if( p->row > maxy ) maxy = p->row;
      if( p->row < miny ) miny = p->row;

    }

    //cout << "(cluster with " << c.vpix.size() << " pixels)" << endl;

    //c.col /= c.vpix.size();
    //c.row /= c.vpix.size();
    c.col /= sumnn; // weighted cluster center
    c.row /= sumnn;

    //c.size = c.vpix.size();
    //c.ncol = maxx-minx+1;
    //c.nrow = maxy-miny+1;
    c.scr = c.vpix.size() + 1024 * ( (maxx-minx+1) + 1024*(maxy-miny+1) ); // compressed size, ncol nrow
    c.mindxy = 999;

    c.vpix.clear(); // save space

    vc.push_back(c); // add cluster to vector

    // look for a new seed = unused pixel:

    while( ( ++seed < pb.size() ) && gone[seed] );

  } // while over seeds

  delete [] gone;

  return vc; // vector of clusters

} // getClus

//------------------------------------------------------------------------------
double GetFWHM( TH1 &h )
{
  int imax = h.GetMaximumBin();
  double ymax = h.GetBinContent( imax );

  double xl = h.GetBinCenter( imax );
  for( int ii = imax; ii > 0; --ii ) { // from peak towards left
    if( h.GetBinContent( ii ) > 0.5*ymax )
      xl = h.GetBinCenter( ii ); // overwritten until found
    else
      break; // against later spikes
  }

  double xr = h.GetBinCenter( imax );
  for( int ii = imax; ii <= h.GetNbinsX(); ++ii ) { // from peak towards right
    if( h.GetBinContent( ii ) > 0.5*ymax )
      xr = h.GetBinCenter( ii ); // overwritten until found
    else
      break; // against later spikes
  }
  double fwhm = xr-xl;
  double dx = h.GetBinWidth( imax );
  if( fwhm > 3.1*dx )
    return fwhm;
  else
    return 3.1*dx;

} // GetFWHM

//------------------------------------------------------------------------------
list < vector <cluster> > oneplane( unsigned ipl, list < vector <pixel> > pxlist )
{
  list < vector < cluster > > clist;

  for( auto ev = pxlist.begin(); ev != pxlist.end(); ++ev ) {

    vector <pixel> pb = *ev;

    vector <cluster> vcl = getClus(pb);

    if( ldbg ) cout << ipl << " clusters " << vcl.size() << endl;

    for( vector<cluster>::iterator cA = vcl.begin(); cA != vcl.end(); ++cA ) {

      // cluster isolation:

      vector<cluster>::iterator cD = cA;
      ++cD;
      for( ; cD != vcl.end(); ++cD ) {
	double dx = cD->col - cA->col;
	double dy = cD->row - cA->row;
	double dxy = sqrt( dx*dx + dy*dy );
	if( dxy < cA->mindxy ) cA->mindxy = dxy;
	if( dxy < cD->mindxy ) cD->mindxy = dxy;
      }

    } // cl A

    clist.push_back(vcl);

  }

  return clist;

} // oneplane

//------------------------------------------------------------------------------
int main( int argc, char* argv[] )
{
  cout << "main " << argv[0] << " called with " << argc << " arguments" << endl;

  if( argc < 2 ) {
    cout << "format: tele run" << endl;
    return 1;
  }

  // run number = last arg

  string runnum( argv[argc-1] );
  int run = atoi( argv[argc-1] );

  cout << "run " << run << endl;
  FileReader * reader;
  if( run < 100 )
    reader = new FileReader( runnum.c_str(), "data/run0000$2R$X" );
  else if( run < 1000 )
    reader = new FileReader( runnum.c_str(), "data/run000$3R$X" );
  else if( run < 10000 )
    reader = new FileReader( runnum.c_str(), "data/run00$4R$X" );
  else if( run < 100000 )
    reader = new FileReader( runnum.c_str(), "data/run0$5R$X" );
  else
    reader = new FileReader( runnum.c_str(), "data/run$6R$X" );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // further arguments:

  int lev = 999222111; // last event

  for( int i = 1; i < argc; ++i ) {

    if( !strcmp( argv[i], "-l" ) )
      lev = atoi( argv[++i] ); // last event

  } // argc

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // runs.dat:

  cout << endl;

  string geoFileName( "geo.dat" );
  double pbeam = 4.8;
  double X0 = 1E-3; // Mimosa or DUT

  ifstream runsFile( "runs.dat" );

  if( runsFile.bad() || ! runsFile.is_open() ) {
    cout << "Error opening runs.dat" << endl;
    return 1;
  }
  else {

    cout << "read runs from runs.dat" << endl;

    string hash( "#" );
    string RUN( "run" );
    string GEO( "geo" );
    string GeV( "GeV" );
    string XNULL( "X0" );
    bool found = 0;

    while( ! runsFile.eof() ) {

      string line;
      getline( runsFile, line );

      if( line.empty() ) continue;

      stringstream tokenizer( line );
      string tag;
      tokenizer >> tag; // leading white space is suppressed
      if( tag.substr(0,1) == hash ) // comments start with #
	continue;

      if( tag == RUN )  {
	int ival;
	tokenizer >> ival;
	if( ival == run ) {
	  found = 1;
	  break; // end file reading
	}
      }

      if( tag == GEO ) {
	tokenizer >> geoFileName;
	continue;
      }

      if( tag == GeV ) {
	tokenizer >> pbeam;
	continue;
      }

      if( tag == XNULL ) {
	tokenizer >> X0;
	continue;
      }

      // anything else on the line and in the file gets ignored

    } // while getline

    if( found )
      cout 
	<< "settings for run " << run << ":" << endl
	<< "  beam " << pbeam << " GeV" << endl
	<< "  geo file " << geoFileName << endl
	<< endl;
    else {
      cout << "run " << run << " not found in runs.dat" << endl;
      return 1;
    }

  } // runsFile

  runsFile.close();

  double ff = 4.8/pbeam;

  const double ang = 0.0024/sqrt(pbeam); // [rad] tritx beam divergence

  const double log10 = log(10);

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // geometry:

  int nx[9]; // x-pixels per plane
  int ny[9]; // y-pixels per plane
  double sizex[9]; // x size per plane
  double sizey[9]; // y size per plane
  double ptchx[9]; // x-pixel size
  double ptchy[9]; // y-pixel size
  double midx[9]; // x mid
  double midy[9]; // y mid

  double zz[9];

  for( int ipl = 0; ipl < 9; ++ipl )
    nx[ipl] = 0; // missing plane flag

  ifstream geoFile( geoFileName );

  if( geoFile.bad() || ! geoFile.is_open() ) {
    cout << "Error opening " << geoFileName << endl;
    return 1;
  }

  cout << "read geometry from " << geoFileName << endl;

  { // open local scope

    string hash{ "#" };
    string plane{ "plane" };
    string type{ "type" };
    string sizexs{ "sizex" };
    string sizeys{ "sizey" };
    string npixelx{ "npixelx" };
    string npixely{ "npixely" };
    string zpos{ "zpos" };

    int ipl = 0;
    string chiptype;

    while( ! geoFile.eof() ) {

      string line;
      getline( geoFile, line );
      cout << line << endl;

      if( line.empty() ) continue;

      stringstream tokenizer( line );
      string tag;
      tokenizer >> tag; // leading white space is suppressed
      if( tag.substr(0,1) == hash ) // comments start with #
	continue;

      if( tag == plane ) {
	tokenizer >> ipl;
	continue;
      }

      if( ipl < 1 || ipl > 6 ) { // Mimosa 1..6
	//cout << "geo wrong plane number " << ipl << endl;
	continue;
      }

      if( tag == type ) {
	tokenizer >> chiptype;
	continue;
      }

      if( tag == sizexs ) {
	double val;
	tokenizer >> val;
	sizex[ipl] = val;
	continue;
      }

      if( tag == sizeys ) {
	double val;
	tokenizer >> val;
	sizey[ipl] = val;
	continue;
      }

      if( tag == npixelx ) {
	int val;
	tokenizer >> val;
	nx[ipl] = val;
	continue;
      }

      if( tag == npixely ) {
	int val;
	tokenizer >> val;
	ny[ipl] = val;
	continue;
      }

      if( tag == zpos ) {
	double val;
	tokenizer >> val;
	zz[ipl] = val;
	continue;
      }

      // anything else on the line and in the file gets ignored

    } // while getline

    for( int ipl = 0; ipl < 9; ++ipl ) {
      if( nx[ipl] == 0 ) continue; // missing plane flag
      ptchx[ipl] = sizex[ipl] / nx[ipl]; // pixel size 21.2/1152=18.4
      ptchy[ipl] = sizey[ipl] / ny[ipl];
      midx[ipl] = 0.5 * sizex[ipl]; // mid plane
      midy[ipl] = 0.5 * sizey[ipl]; // mid plane
    }

  } // geo scope

  geoFile.close();

  double scat = 0.0136 * sqrt(X0) / pbeam * ( 1 + 0.038*log(X0) ); // [rad] DUT

  cout << endl << "DUT scat " << scat*1E3 << " mrad" << endl;

  double DUTz = 0.5 * ( zz[3] + zz[4] );
  double dzd = zz[6] - DUTz; // driplet lever arm
  double trng = scat * dzd; // DUT scattering
  cout << "fiducial " << 3*trng << " mm " << endl;

  double X0m = 0.9E-3; // Mimosa and air
  double scatm = 0.0136 * sqrt(X0m) / pbeam * ( 1 + 0.038*log(X0m) ); // [rad] Mimosa scattering

  cout << endl << "Mim scat " << scatm*1E3 << " mrad" << endl;

  // add resolution:

  double dzt = zz[2] - zz[1];
  scatm += 0.005 / dzt;

  cout << endl << "kink width " << scatm*1E3 << " mrad" << endl;

  double tricut = 3 * scatm; // [rad]

  cout << endl << "kinkCut " << tricut*1E3 << " mrad" << endl;

  double rmsm = 0.005 + scatm * dzt; // [mm] Mimosa

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // Mimosa hot pixels:

  ostringstream hotFileName; // output string stream

  hotFileName << "hot_" << run << ".dat";

  ifstream ihotFile( hotFileName.str() );

  set <int> hotset[9];

  cout << endl;
  if( ihotFile.bad() || ! ihotFile.is_open() ) {
    cout << "no " << hotFileName.str() << ", will be created" << endl;
  }
  // can there be instructions between if and else ? no!
  else {

    cout << "read hot pixel list from " << hotFileName.str() << endl;

    string hash{ "#" };
    string plane{ "plane" };
    string pix{ "pix" };

    int ipl = 0;

    while( ! ihotFile.eof() ) {

      string line;
      getline( ihotFile, line );
      //cout << line << endl;

      if( line.empty() ) continue;

      stringstream tokenizer( line );
      string tag;
      tokenizer >> tag; // leading white space is suppressed
      if( tag.substr(0,1) == hash ) // comments start with #
	continue;

      if( tag == plane )
	tokenizer >> ipl;

      if( ipl < 1 || ipl > 6 ) {
	cout << "hot pixel wrong plane number " << ipl << endl;
	continue;
      }

      if( tag == pix ) {
	int ix, iy;
	tokenizer >> ix;
	tokenizer >> iy;
	int ipx = ix*ny[ipl]+iy;
	hotset[ipl].insert(ipx);
      }

    } // while getline

  } // hotFile

  ihotFile.close();

  for( int ipl = 0; ipl < 9; ++ipl )
    cout << ipl << ": hot " << hotset[ipl].size() << endl;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // telescope alignment:

  int aligniteration = 0;
  double alignx[9];
  double aligny[9];
  double alignz[9];
  double rotx[9];
  double roty[9];

  for( int ipl = 0; ipl < 9; ++ipl ) {

    alignx[ipl] = 0.000; // [mm] same sign as dxAB
    aligny[ipl] = 0.000; // [mm] same sign as dy
    alignz[ipl] = 0.000; // [mm]
    rotx[ipl] = 0.0000; // [rad] rot, same     sign dxvsy
    roty[ipl] = 0.0000; // [rad] rot, opposite sign dyvsx

  }

  ostringstream alignFileName; // output string stream

  alignFileName << "align_" << run << ".dat";

  ifstream ialignFile( alignFileName.str() );

  cout << endl;
  if( ialignFile.bad() || ! ialignFile.is_open() ) {
    cout << "no " << alignFileName.str() << ", will bootstrap" << endl;
    cout << endl;
  }
  else {

    cout << "read alignment from " << alignFileName.str() << endl;

    string hash{ "#" };
    string iteration{ "iteration" };
    string plane{ "plane" };
    string shiftx{ "shiftx" };
    string shifty{ "shifty" };
    string shiftz{ "shiftz" };
    string rotxvsy{ "rotxvsy" };
    string rotyvsx{ "rotyvsx" };

    int ipl = 0;

    while( ! ialignFile.eof() ) {

      string line;
      getline( ialignFile, line );
      cout << line << endl;

      if( line.empty() ) continue;

      stringstream tokenizer( line );
      string tag;
      tokenizer >> tag; // leading white space is suppressed
      if( tag.substr(0,1) == hash ) // comments start with #
	continue;

      if( tag == iteration ) {
	tokenizer >> aligniteration;
	continue;
      }

      if( tag == plane )
	tokenizer >> ipl;

      if( ipl < 1 || ipl > 6 ) {
	cout << "align wrong plane number " << ipl << endl;
	continue;
      }

      double val;
      tokenizer >> val;
      if(      tag == shiftx )
	alignx[ipl] = val;
      else if( tag == shifty )
	aligny[ipl] = val;
      else if( tag == shiftz )
	alignz[ipl] = val;
      else if( tag == rotxvsy )
	rotx[ipl] = val;
      else if( tag == rotyvsx )
	roty[ipl] = val;

      // anything else on the line and in the file gets ignored

    } // while getline

  } // alignFile

  ialignFile.close();

  if( aligniteration == 0 ) ff *= 3; // wider binning

  double z2 = zz[2] + alignz[2];
  double z3 = zz[3] + alignz[3];
  double dz32 = z3 - z2;
  double zm23 = ( z3 + z2 ) * 0.5;

  double z4 = zz[4] + alignz[4];
  double z5 = zz[5] + alignz[5];
  double dz54 = z5 - z4;
  double zm45 = ( z5 + z4 ) * 0.5;

  //double dz43 = z4 - z3;
  //double sixcut = 3 * ( 0.005 + 0.5*scat*dz43 ); // [mm]
  double sixcut = 0.025; // [mm] zt z intersect
  if( run == 36941 )
    sixcut = 0.165; // lots of air
  if( run == 36942 )
    sixcut = 0.300; // lots of air
  if( run == 36943 )
    sixcut = 0.400; // lots of air
  if( run == 36944 )
    sixcut = 0.400; // lots of air
  if( run == 36945 )
    sixcut = 0.600; // lots of air

  cout << endl << "sixcut " << sixcut*1E3 << " um" << endl;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // (re-)create root file:

  ostringstream rootFileName; // output string stream

  rootFileName << "tele" << run << ".root";

  TFile* histoFile = new TFile( rootFileName.str(  ).c_str(  ), "RECREATE" );

  // book histos:

  TH1I hdtus = TH1I( "dtus", "time between events;time between events [us];events", 100, 0, 1000 );
  TH1I hdtms = TH1I( "dtms", "time between events;time between events [ms];events", 100, 0, 1000 );

  TH1I t1Histo = TH1I( "t1", "event time;event time [s];events / 10 ms", 100, 0, 1 );
  TH1I t2Histo = TH1I( "t2", "event time;event time [s];events / s", 300, 0, 300 );
  TH1I t3Histo = TH1I( "t3", "event time;event time [s];events / 10 s", 150, 0, 1500 );
  TH1I t4Histo = TH1I( "t4", "event time;event time [s];events / 10 s", 600, 0, 6000 );
  TH1I t5Histo = TH1I( "t5", "event time;event time [s];events / min", 1000, 0, 60000 );

  TH1I hnpx0[9];
  TH1I hpivot[9];
  TH1I hcol0[9];
  TH1I hrow0[9];
  TH2I * hmap0[9];

  TH1I hcol[9];
  TH1I hrow[9];
  TH1I hbool[9];
  TH1I hnpx[9];

  TH1I hncl[9];
  TH1I hccol[9];
  TH1I hcrow[9];
  TH1I hnpix[9];
  TH1I hncol[9];
  TH1I hnrow[9];
  TH1I hmindxy[9];

  TH2I * hxx[9];
  TH1I hdx[9];
  TH1I hdy[9];
  TProfile dxvsx[9];
  TProfile dxvsy[9];
  TProfile dyvsx[9];
  TProfile dyvsy[9];

  TH1I hdx4[9];
  TH1I hdy4[9];
  TProfile dx4vsy[9];
  TProfile dy4vsx[9];

  TProfile effvsx[9];
  TProfile effvsy[9];
  TProfile effvsxm[9];
  TProfile effvsym[9];
  TProfile2D * effvsxmym[9];

  for( int ipl = 1; ipl <= 6; ++ipl ) {

    hpivot[ipl] = TH1I( Form( "pivot%i", ipl ),
			Form( "plane %i pivot;pivot;plane %i events", ipl, ipl ),
			800, 0, 8000 );
    hnpx0[ipl] = TH1I( Form( "allnpx%i", ipl ),
		       Form( "plane %i all pixels per event;all pixels / event;plane %i events", ipl, ipl ),
		       200, 0, 200 );
    hcol0[ipl] = TH1I( Form( "allcol%i", ipl ),
		       Form( "plane %i all pix col;col;all plane %i pixels", ipl, ipl ), 
		       nx[ipl]/4, 0, nx[ipl] );
    hrow0[ipl] = TH1I( Form( "allrow%i", ipl ),
		       Form( "plane %i all pix row;row;all plane %i pixels", ipl, ipl ),
		       ny[ipl]/2, 0, ny[ipl] );
    hmap0[ipl] = new TH2I( Form( "map%i", ipl ),
			   Form( "plane %i map all;col;row;all plane %i pixels", ipl, ipl ),
			   nx[ipl]/4, 0, nx[ipl], ny[ipl]/2, 0, ny[ipl] );

    hnpx[ipl] = TH1I( Form( "npx%i", ipl ),
		      Form( "plane %i cool pixel per event;cool pixels / event;plane %i events", ipl, ipl ),
		      200, 0, 200 );
    hcol[ipl] = TH1I( Form( "col%i", ipl ),
		      Form( "plane %i cool pix col;col;cool plane %i pixels", ipl, ipl ), 
		      nx[ipl]/4, 0, nx[ipl] );
    hrow[ipl] = TH1I( Form( "row%i", ipl ),
		      Form( "plane %i cool pix row;row;cool plane %i pixels", ipl, ipl ),
		      ny[ipl]/2, 0, ny[ipl] );
    hbool[ipl] = TH1I( Form( "bool%i", ipl ),
		       Form( "plane %i before/after pivot;before/after pivot;cool plane %i pixels", ipl, ipl ),
		       2, -0.4, 1.5 );

    hncl[ipl] = TH1I( Form( "ncl%i", ipl ),
		      Form( "plane %i cluster per event;cluster;plane %i events", ipl, ipl ),
		      51, -0.5, 50.5 );
    hccol[ipl] = TH1I( Form( "ccol%i", ipl ),
		       Form( "plane %i cluster col;col;%i clusters", ipl, ipl ), 
		       nx[ipl]/4, 0, nx[ipl] );
    hcrow[ipl] = TH1I( Form( "crow%i", ipl ),
		       Form( "plane %i cluster row;row;%i clusters", ipl, ipl ),
		       ny[ipl]/2, 0, ny[ipl] );
    hnpix[ipl] = TH1I( Form( "npix%i", ipl ),
		       Form( "plane %i cluster size;pixels/cluster;%i clusters", ipl, ipl ),
		       80, 0.5, 80.5 );
    hncol[ipl] = TH1I( Form( "ncol%i", ipl ), 
		       Form( "plane %i cluster size x;columns/cluster;%i clusters", ipl, ipl ),
		       50, 0.5, 50.5 );
    hnrow[ipl] = TH1I( Form( "nrow%i", ipl ),
		       Form( "plane %i cluster size y;rows/cluster;%i clusters", ipl, ipl ),
		       50, 0.5, 50.5 );

    hmindxy[ipl] = TH1I( Form( "mindxy%i", ipl ),
			 Form( "plane %i cluster isolation;distance to next cluster [pixels];%i clusters",
			       ipl, ipl ),
			 100, 0, 10 );

    hxx[ipl] = new TH2I( Form( "xx%im", ipl ),
			 Form( "plane %i-m x-x;mid  plane x [mm];plane %i x [mm];cluster pairs", ipl, ipl ),
			 nx[ipl]/4, -midx[ipl], midx[ipl], nx[ipl]/4, -midx[ipl], midx[ipl] );

    hdx[ipl] = TH1I( Form( "dx%im", ipl ),
		     Form( "plane %i-m dx;%i-m dx [mm];cluster pairs", ipl, ipl ),
		     400, -2, 2 );
    hdy[ipl] = TH1I( Form( "dy%im", ipl ),
		     Form( "plane %i-m dy;%i-m dy [mm];cluster pairs", ipl, ipl ),
		     400, -2, 2 );

    dxvsx[ipl] = TProfile( Form( "dx%imvsx", ipl ),
			   Form( "plane %i-m dx vs x;x [mm];<%i-m dx> [mm]", ipl, ipl ),
			   100, -midx[ipl], midx[ipl], -0.1, 0.1 );
    dxvsy[ipl] = TProfile( Form( "dx%imvsy", ipl ),
			   Form( "plane %i-m dx vs y;y [mm];<%i-m dx> [mm]", ipl, ipl ),
			   100, -midy[ipl], midy[ipl], -0.1, 0.1 );
    dyvsx[ipl] = TProfile( Form( "dy%imvsx", ipl ),
			   Form( "plane %i-m dy vs x;x [mm];<%i-m dy> [mm]", ipl, ipl ),
			   100, -midx[ipl], midx[ipl], -0.1, 0.1 );
    dyvsy[ipl] = TProfile( Form( "dy%imvsy", ipl ),
			   Form( "plane %i-m dy vs y;y [mm];<%i-m dy> [mm]", ipl, ipl ),
			   100, -midy[ipl], midy[ipl], -0.1, 0.1 );

    hdx4[ipl] = TH1I( Form( "dx4%i", ipl ),
		      Form( "plane %i dx4;plane %i dx4 [#mum];track-cluster pairs", ipl, ipl ),
		      200, -100, 100 );
    hdy4[ipl] = TH1I( Form( "dy4%i", ipl ),
		      Form( "plane %i dy4;plane %i dy4 [#mum];track-cluster pairs", ipl, ipl ),
		      200, -100, 100 );
    dx4vsy[ipl] = TProfile( Form( "dx4%ivsy", ipl ),
			    Form( "plane %i dx4 vs y;y [mm];plane %i <dx4> [#mum]", ipl, ipl ),
			    100, -midy[ipl], midy[ipl], -100, 100 );
    dy4vsx[ipl] = TProfile( Form( "dy4%ivsx", ipl ),
			    Form( "plane %i dy4 vs x;x [mm];plane %i <dy4> [#mum]", ipl, ipl ),
			    100, -midx[ipl], midx[ipl], -100, 100 );

    effvsx[ipl] = TProfile( Form( "eff%ivsx", ipl ),
			    Form( "plane %i efficiency;x [mm];plane %i efficiency", ipl, ipl ),
			    100, -midx[ipl], midx[ipl], -1, 2 );
    effvsy[ipl] = TProfile( Form( "eff%ivsy", ipl ),
			    Form( "plane %i efficiency;y [mm];plane %i efficiency", ipl, ipl ),
			    100, -midy[ipl], midy[ipl], -1, 2 );
    effvsxm[ipl] = TProfile( Form( "eff%ivsxm", ipl ),
			     Form( "plane %i efficiency;x mod 36.8 [#mum];plane %i efficiency", ipl, ipl ),
			     74, 0, 2*ptchx[ipl]*1E3, -1, 2 );
    effvsym[ipl] = TProfile( Form( "eff%ivsym", ipl ),
			     Form( "plane %i efficiency;y mod 36.8 [#mum];plane %i efficiency", ipl, ipl ),
			     74, 0, 2*ptchy[ipl]*1E3, -1, 2 );
    effvsxmym[ipl] = new
      TProfile2D( Form( "eff%ivsxmym", ipl ),
		  Form( "plane %i efficiency;x mod 36.8 [#mum];y mod 36.8 [#mum];plane %i efficiency", ipl, ipl ),
		  74, 0, 2*ptchx[ipl]*1E3, 74, 0, 2*ptchy[ipl]*1E3, -1, 2 );

  } // planes

  TH1I hdxCA[2];
  TH1I hdyCA[2];
  TProfile dxCAvsx[2];
  TProfile dyCAvsy[2];

  TH1I htridx[2];
  TH1I htridy[2];

  TH1I htrikx[2];
  TH1I htrikxc[2];
  TH1I htrikxci[2];
  TH1I htriky[2];
  TH1I htrikyc[2];
  TH1I htrikyci[2];
  TProfile trimadkxvsx[2];
  TProfile trimadkxvstx[2];

  TH1I htridxc[2];
  TH1I htridxci[2];
  TH1I htridxct[2];
  TH1I htridxctg[2];

  TH1I htrixm[2];
  TH1I htrixm1[2];
  TH1I htrixm2[2];
  TH1I htrixm3[2];
  TH1I htrixm4[2];

  TH1I htridxc1[2];
  TH1I htridxc2[2];
  TH1I htridxc3[2];
  TH1I htridxc4[2];
  TH1I htridxc5[2];
  TH1I htridxc6[2];

  TH1I htridxc21[2];
  TH1I htridxc22[2];
  TH1I htridxc23[2];
  TH1I htridxc24[2];
  TH1I htridxc25[2];
  TH1I htridxc223[2];
  TH1I htridxc224[2];
  TH1I htridxctgg[2];
  TH1I htridxc111[2];

  TProfile tridxvsx[2];
  TProfile tridxvsy[2];
  TProfile tridxvsxm[2];
  TProfile trimadxvsxm[2];
  TProfile2D * trimadxvsxmym[2];
  TProfile tridxvstx[2];
  TProfile trimadxvstx[2];

  TH1I htridyc[2];
  TH1I htridyci[2];
  TH1I htridyct[2];
  TH1I htridyc1[2];
  TH1I htridyc2[2];
  TH1I htridyc3[2];
  TH1I htridyc4[2];
  TH1I htridyc5[2];
  TH1I htridyc6[2];
  TH1I htridycg[2];

  TProfile tridyvsx[2];
  TProfile tridyvsy[2];
  TProfile tridyvsym[2];
  TProfile trimadyvsym[2];
  TProfile2D * trimadyvsxmym[2];
  TProfile tridyvsty[2];
  TProfile trimadyvsty[2];

  TProfile tridxCvsx[2];
  TProfile tridxCvsy[2];
  TProfile tridyCvsx[2];
  TProfile tridyCvsy[2];
  TProfile tridxCvsax[2];
  TProfile tridyCvsay[2];

  TH1I htrix[2];
  TH1I htriy[2];
  TH2I * htrixy[2];
  TH1I htritx[2];
  TH1I htrity[2];

  TH1I htrincol[2];
  TH1I htrinrow[2];
  TH1I htrinpix[2];
  TH2D * hnpx1map[2];
  TH2D * hnpx2map[2];
  TH2D * hnpx3map[2];
  TH2D * hnpx4map[2];
  TProfile trincolvsxm[2];
  TProfile trinrowvsym[2];
  TProfile2D * trinpixvsxmym[2];
  TProfile2D * trinpixgvsxmym[2];

  for( int itd = 0; itd < 2; ++itd ) {

    string tri{"tri"};
    string dri{"dri"};
    string tds{tri};
    if( itd ) 
      tds = dri;
    int ipl = 2+3*itd;

    hdxCA[itd] = TH1I( Form( "%sdxCA", tds.c_str() ),
		       Form( "%splet dx CA;%splet dx CA [mm];C-A pairs",
			     tds.c_str(), tds.c_str() ),
		       100, -1, 1 );
    hdyCA[itd] = TH1I( Form( "%sdyCA", tds.c_str() ),
		       Form( "%splet dy CA;%splet dy CA [mm];C-A pairs",
			     tds.c_str(), tds.c_str() ),
		       100, -1, 1 );

    dxCAvsx[itd] = TProfile( Form( "%sdxCAvsx", tds.c_str() ),
			     Form( "%splet dx CA vs x;%splet xC [mm];<%splets #Deltax CA> [mm]",
				   tds.c_str(), tds.c_str(), tds.c_str() ),
			     100, -midx[ipl], midx[ipl], -1, 1 );
    dyCAvsy[itd] = TProfile( Form( "%sdyCAvsy", tds.c_str() ),
			     Form( "%splet dy CA vs y;%splet yC [mm];<%splets #Deltay CA> [mm]",
				   tds.c_str(), tds.c_str(), tds.c_str() ),
			     100, -midy[ipl], midy[ipl], -1, 1 );

    htridx[itd] = TH1I( Form( "%sdx", tds.c_str() ),
			Form( "%splet dx;%splet dx [#mum];%splets",
			      tds.c_str(), tds.c_str(), tds.c_str() ),
			200, -10E3*rmsm, 10E3*rmsm );

    htridy[itd] = TH1I( Form( "%sdy", tds.c_str() ),
			Form( "%splet dy;%splet dy [#mum];%splets",
			      tds.c_str(), tds.c_str(), tds.c_str() ),
			200, -10E3*rmsm, 10E3*rmsm );

    htrikx[itd] =  TH1I( Form( "%skx", tds.c_str() ),
			 Form( "%splet kink x;kink x [mrad];%splets",
			       tds.c_str(), tds.c_str() ),
			 200, -10E3*scatm, 10E3*scatm );
    htriky[itd] =  TH1I( Form( "%sky", tds.c_str() ),
			 Form( "%splet kink y;kink y [mrad];%splets",
			       tds.c_str(), tds.c_str() ),
			 200, -10E3*scatm, 10E3*scatm );
    htrikxc[itd] =  TH1I( Form( "%skxc", tds.c_str() ),
			  Form( "%splet kink x;kink x [mrad];%splets",
				tds.c_str(), tds.c_str() ),
			  200, -10E3*scatm, 10E3*scatm );
    htrikxci[itd] =  TH1I( Form( "%skxci", tds.c_str() ),
			   Form( "%splet kink x;kink x [mrad];isolated %splets",
				 tds.c_str(), tds.c_str() ),
			   200, -10E3*scatm, 10E3*scatm );
    htrikyc[itd] =  TH1I( Form( "%skyc", tds.c_str() ),
			  Form( "%splet kink y;kink y [mrad];%splets",
				tds.c_str(), tds.c_str() ),
			  200, -10E3*scatm, 10E3*scatm );
    htrikyci[itd] =  TH1I( Form( "%skyci", tds.c_str() ),
			   Form( "%splet kink y;kink y [mrad];isolated %splets",
				 tds.c_str(), tds.c_str() ),
			   200, -10E3*scatm, 10E3*scatm );

    trimadkxvsx[itd] = TProfile( Form( "%smadkxvsx", tds.c_str() ),
				 Form( "%splet MAD kx vs x;x [mm];%splet MAD(kx) [mrad]",
				       tds.c_str(), tds.c_str() ),
				 110, -11, 11, 0, 5E3*scatm );
    trimadkxvstx[itd] = TProfile( Form( "%smadkxvstx", tds.c_str() ),
				  Form( "%splet MAD kx vs x;slope x [mrad];%splet MAD(kx) [mrad]",
					tds.c_str(), tds.c_str() ),
				  80, -3E3*(ang+itd*scat), 3E3*(ang+itd*scat), 0, 5E3*scatm );

    htridxc[itd] = TH1I( Form( "%sdxc", tds.c_str() ),
			 Form( "%splet dx;%splet dx [#mum];%splets",
			       tds.c_str(), tds.c_str(), tds.c_str() ),
			 200, -10E3*rmsm, 10E3*rmsm );
    htridxci[itd] = TH1I( Form( "%sdxci", tds.c_str() ),
			  Form( "isolated %splet dx;%splet dx [#mum];isolated %splets",
				tds.c_str(), tds.c_str(), tds.c_str() ),
			  200, -10E3*rmsm, 10E3*rmsm );

    htridxct[itd] = TH1I( Form( "%sdxct", tds.c_str() ),
			  Form( "%splet dx;%splet dx [#mum];%splets",
				tds.c_str(), tds.c_str(), tds.c_str() ),
			  200, -10E3*rmsm, 10E3*rmsm );
    htridxctg[itd] = TH1I( Form( "%sdxctg", tds.c_str() ),
			   Form( "good %splet dx;%splet dx [#mum];good col %splets",
				 tds.c_str(), tds.c_str(), tds.c_str() ),
			   200, -10E3*rmsm, 10E3*rmsm );
    htridxctgg[itd] = TH1I( Form( "%sdxctgg", tds.c_str() ),
			    Form( "%splet dx good col;%splet dx [#mum];good-col %splets",
				  tds.c_str(), tds.c_str(), tds.c_str() ),
			    200, -10E3*rmsm, 10E3*rmsm );

    htridxc1[itd] = TH1I( Form( "%sdxc1", tds.c_str() ),
			  Form( "%splet dx 1-col;%splet dx [#mum];1-col %splets",
				tds.c_str(), tds.c_str(), tds.c_str() ),
			  200, -10E3*rmsm, 10E3*rmsm );
    htridxc2[itd] = TH1I( Form( "%sdxc2", tds.c_str() ),
			  Form( "%splet dx 2-col;%splet dx [#mum];2-col %splets",
				tds.c_str(), tds.c_str(), tds.c_str() ),
			  200, -10E3*rmsm, 10E3*rmsm );
    htridxc3[itd] = TH1I( Form( "%sdxc3", tds.c_str() ),
			  Form( "%splet dx 3-col;%splet dx [#mum];3-col %splets",
				tds.c_str(), tds.c_str(), tds.c_str() ),
			  200, -10E3*rmsm, 10E3*rmsm );
    htridxc4[itd] = TH1I( Form( "%sdxc4", tds.c_str() ),
			  Form( "%splet dx 4-col;%splet dx [#mum];4-col %splets",
				tds.c_str(), tds.c_str(), tds.c_str() ),
			  200, -10E3*rmsm, 10E3*rmsm );
    htridxc5[itd] = TH1I( Form( "%sdxc5", tds.c_str() ),
			  Form( "%splet dx 5-col;%splet dx [#mum];5-col %splets",
				tds.c_str(), tds.c_str(), tds.c_str() ),
			  200, -10E3*rmsm, 10E3*rmsm );
    htridxc6[itd] = TH1I( Form( "%sdxc6", tds.c_str() ),
			  Form( "%splet dx 6-col;%splet dx [#mum];6-col %splets",
				tds.c_str(), tds.c_str(), tds.c_str() ),
			  200, -10E3*rmsm, 10E3*rmsm );
    htridxc21[itd] = TH1I( Form( "%sdxc21", tds.c_str() ),
			   Form( "%splet dx 2-col 1-row;%splet dx [#mum];2-col 1-row %splets",
				 tds.c_str(), tds.c_str(), tds.c_str() ),
			   200, -10E3*rmsm, 10E3*rmsm );
    htridxc22[itd] = TH1I( Form( "%sdxc22", tds.c_str() ),
			   Form( "%splet dx 2-col 2-row;%splet dx [#mum];2-col 2-row %splets",
				 tds.c_str(), tds.c_str(), tds.c_str() ),
			   200, -10E3*rmsm, 10E3*rmsm );
    htridxc23[itd] = TH1I( Form( "%sdxc23", tds.c_str() ),
			   Form( "%splet dx 2-col 3-row;%splet dx [#mum];2-col 3-row %splets",
				 tds.c_str(), tds.c_str(), tds.c_str() ),
			   200, -10E3*rmsm, 10E3*rmsm );
    htridxc24[itd] = TH1I( Form( "%sdxc24", tds.c_str() ),
			   Form( "%splet dx 2-col 4-row;%splet dx [#mum];2-col 4-row %splets",
				 tds.c_str(), tds.c_str(), tds.c_str() ),
			   200, -10E3*rmsm, 10E3*rmsm );
    htridxc25[itd] = TH1I( Form( "%sdxc25", tds.c_str() ),
			   Form( "%splet dx 2-col 5-row;%splet dx [#mum];2-col 5-row %splets",
				 tds.c_str(), tds.c_str(), tds.c_str() ),
			   200, -10E3*rmsm, 10E3*rmsm );
    htridxc223[itd] = TH1I( Form( "%sdxc223", tds.c_str() ),
			    Form( "%splet dx 2-col 2-row 3-pix;%splet dx [#mum];2-col 2-row 3-pix %splets",
				  tds.c_str(), tds.c_str(), tds.c_str() ),
			    200, -10E3*rmsm, 10E3*rmsm );
    htridxc224[itd] = TH1I( Form( "%sdxc224", tds.c_str() ),
			    Form( "%splet dx 2-col 2-row 4-pix;%splet dx [#mum];2-col 2-row 4-pix %splets",
				  tds.c_str(), tds.c_str(), tds.c_str() ),
			    200, -10E3*rmsm, 10E3*rmsm );
    htridxc111[itd] = TH1I( Form( "%sdxc111", tds.c_str() ),
			    Form( "%splet dx 111-col;111-col %splet dx [#mum];111-col %splets",
				  tds.c_str(), tds.c_str(), tds.c_str() ),
			    200, -10E3*rmsm, 10E3*rmsm );

    tridxvsx[itd] = TProfile( Form( "%sdxvsx", tds.c_str() ),
			      Form( "%splet dx vs x;%splet xB [mm];<%splets #Deltax> [#mum]",
				    tds.c_str(), tds.c_str(), tds.c_str() ),
			      100, -midx[ipl], midx[ipl], -50, 50 );
    tridxvsy[itd] = TProfile( Form( "%sdxvsy", tds.c_str() ),
			      Form( "%splet dx vs y;%splet yB [mm];<%splets #Deltax> [#mum]",
				    tds.c_str(), tds.c_str(), tds.c_str() ),
			      100, -midy[ipl], midy[ipl], -50, 50 );
    tridxvsxm[itd] = TProfile( Form( "%sdxvsxm", tds.c_str() ),
			       Form( "%splet dx vs xmod;%splet xB mod 36.8 [#mum];<%splets #Deltax> [#mum]",
				     tds.c_str(), tds.c_str(), tds.c_str() ),
			       74, 0, 2*ptchx[ipl]*1E3, -50, 50 );
    trimadxvsxm[itd] =
      TProfile( Form( "%smadxvsxm", tds.c_str() ),
		Form( "%splet MAD(#Deltax) vs xmod;%splet xB mod 36.8;%splet MAD(#Deltax) [#mum]",
		      tds.c_str(), tds.c_str(), tds.c_str() ),
		74, 0, 2*ptchx[ipl]*1E3, 0, 50 );

    trimadxvsxmym[itd] = new
      TProfile2D( Form( "%smadxvsxmym", tds.c_str() ),
		  Form( "%splet MAD(#Deltax) vs xmod;%splet xB mod 36.8;%splet yB mod 36.8;%splet MAD(#Deltax) [#mum]",
			tds.c_str(), tds.c_str(), tds.c_str(), tds.c_str() ),
		  74, 0, 2*ptchx[ipl]*1E3, 74, 0, 2*ptchy[ipl]*1E3, 0, 50 );

    tridxvstx[itd] = TProfile( Form( "%sdxvstx", tds.c_str() ),
			       Form( "%splet dx vs tx;%splet slope x [mrad];<%splets #Deltax> [#mum]",
				     tds.c_str(), tds.c_str(), tds.c_str() ),
			       80, -3E3*(ang+itd*scat), 3E3*(ang+itd*scat), -50, 50 );
    trimadxvstx[itd] =
      TProfile( Form( "%smadxvstx", tds.c_str() ),
		Form( "%splet MAD(#Deltax) vs #theta_{x};%splet #theta_{x} [mrad];%splet MAD(#Deltax) [#mum]",
		      tds.c_str(), tds.c_str(), tds.c_str() ),
		100, -3E3*(ang+itd*scat), 3E3*(ang+itd*scat), 0, 50 );

    htrixm[itd] = TH1I( Form( "%sxm", tds.c_str() ),
			Form( "%splet xmod;%splet x mod 18.4 [#mum];%splets",
			      tds.c_str(), tds.c_str(), tds.c_str() ),
			74, 0, 1*ptchx[ipl]*1E3 );
    htrixm1[itd] = TH1I( Form( "%sxm1", tds.c_str() ),
			 Form( "%splet xmod 1-col;%splet x mod 18.4 [#mum];1-col %splets",
			       tds.c_str(), tds.c_str(), tds.c_str() ),
			 74, 0, 1*ptchx[ipl]*1E3 );
    htrixm2[itd] = TH1I( Form( "%sxm2", tds.c_str() ),
			 Form( "%splet xmod 2-col;%splet x mod 18.4 [#mum];2-col %splets",
			       tds.c_str(), tds.c_str(), tds.c_str() ),
			 74, 0, 1*ptchx[ipl]*1E3 );
    htrixm3[itd] = TH1I( Form( "%sxm3", tds.c_str() ),
			 Form( "%splet xmod 3-col;%splet x mod 18.4 [#mum];3-col %splets",
			       tds.c_str(), tds.c_str(), tds.c_str() ),
			 74, 0, 1*ptchx[ipl]*1E3 );
    htrixm4[itd] = TH1I( Form( "%sxm4", tds.c_str() ),
			 Form( "%splet xmod 4-col;%splet x mod 18.4 [#mum];4-col %splets",
			       tds.c_str(), tds.c_str(), tds.c_str() ),
			 74, 0, 1*ptchx[ipl]*1E3 );

    htridyc[itd] = TH1I( Form( "%sdyc", tds.c_str() ),
			 Form( "%splet dy;%splet dy [#mum];%splets",
			       tds.c_str(), tds.c_str(), tds.c_str() ),
			 200, -10E3*rmsm, 10E3*rmsm );
    htridyci[itd] = TH1I( Form( "%sdyci", tds.c_str() ),
			  Form( "isolated %splet dy;%splet dy [#mum];isolated %splets",
				tds.c_str(), tds.c_str(), tds.c_str() ),
			  200, -10E3*rmsm, 10E3*rmsm );
    htridyct[itd] = TH1I( Form( "%sdyct", tds.c_str() ),
			  Form( "%splet dy;%splet dy [#mum];%splets",
				tds.c_str(), tds.c_str(), tds.c_str() ),
			  200, -10E3*rmsm, 10E3*rmsm );

    htridyc1[itd] = TH1I( Form( "%sdyc1", tds.c_str() ),
			  Form( "%splet dy 1-row;%splet dy [#mum];1-row %splets",
				tds.c_str(), tds.c_str(), tds.c_str() ),
			  200, -10E3*rmsm, 10E3*rmsm );
    htridyc2[itd] = TH1I( Form( "%sdyc2", tds.c_str() ),
			  Form( "%splet dy 2-row;%splet dy [#mum];2-row %splets",
				tds.c_str(), tds.c_str(), tds.c_str() ),
			  200, -10E3*rmsm, 10E3*rmsm );
    htridyc3[itd] = TH1I( Form( "%sdyc3", tds.c_str() ),
			  Form( "%splet dy 3-row;%splet dy [#mum];3-row %splets",
				tds.c_str(), tds.c_str(), tds.c_str() ),
			  200, -10E3*rmsm, 10E3*rmsm );
    htridyc4[itd] = TH1I( Form( "%sdyc4", tds.c_str() ),
			  Form( "%splet dy 4-row;%splet dy [#mum];4-row %splets",
				tds.c_str(), tds.c_str(), tds.c_str() ),
			  200, -10E3*rmsm, 10E3*rmsm );
    htridyc5[itd] = TH1I( Form( "%sdyc5", tds.c_str() ),
			  Form( "%splet dy 5-row;%splet dy [#mum];5-row %splets",
				tds.c_str(), tds.c_str(), tds.c_str() ),
			  200, -10E3*rmsm, 10E3*rmsm );
    htridyc6[itd] = TH1I( Form( "%sdyc6", tds.c_str() ),
			  Form( "%splet dy 6-row;%splet dy [#mum];6-row %splets",
				tds.c_str(), tds.c_str(), tds.c_str() ),
			  200, -10E3*rmsm, 10E3*rmsm );
    htridycg[itd] = TH1I( Form( "%sdycg", tds.c_str() ),
			  Form( "%splet dy good row;%splet dy [#mum];good-row %splets",
				tds.c_str(), tds.c_str(), tds.c_str() ),
			  200, -10E3*rmsm, 10E3*rmsm );

    tridyvsx[itd] = TProfile( Form( "%sdyvsx", tds.c_str() ),
			      Form( "%splet dy vs x;%splet xB [mm];<%splets #Deltay> [#mum]",
				    tds.c_str(), tds.c_str(), tds.c_str() ),
			      100, -midx[ipl], midx[ipl], -50, 50 );
    tridyvsy[itd] = TProfile( Form( "%sdyvsy", tds.c_str() ),
			      Form( "%splet dy vs y;%splet yB [mm];<%splets #Deltay> [#mum]",
				    tds.c_str(), tds.c_str(), tds.c_str() ),
			      100, -midy[ipl], midy[ipl], -50, 50 );
    tridyvsym[itd] = TProfile( Form( "%sdyvsym", tds.c_str() ),
			       Form( "%splet dy vs ymod;%splet yB mod 36.8 [#mum];<%splets #Deltay> [#mum]",
				     tds.c_str(), tds.c_str(), tds.c_str() ),
			       74, 0, 2*ptchy[ipl]*1E3, -50, 50 );
    trimadyvsym[itd] =
      TProfile( Form( "%smadyvsym", tds.c_str() ),
		Form( "%splet MAD(#Deltay) vs ymod;%splet yB mod 36.8;%splet MAD(#Deltay) [#mum]",
		      tds.c_str(), tds.c_str(), tds.c_str() ),
		74, 0, 2*ptchy[ipl]*1E3, 0, 50 );

    trimadyvsxmym[itd] = new
      TProfile2D( Form( "%smadyvsxmym", tds.c_str() ),
		  Form( "%splet MAD(#Deltay) vs xmod ymod;%splet xB mod 36.8;%splet yB mod 36.8;%splet MAD(#Deltay) [#mum]",
			tds.c_str(), tds.c_str(), tds.c_str(), tds.c_str() ),
		  74, 0, 2*ptchx[ipl]*1E3, 74, 0, 2*ptchy[ipl]*1E3, 0, 50 );

    tridyvsty[itd] = TProfile( Form( "%sdyvsty", tds.c_str() ),
			       Form( "%splet dy vs ty;%splet slope y [mrad];<%splets #Deltay> [#mum]",
				     tds.c_str(), tds.c_str(), tds.c_str() ),
			       80, -3E3*(ang+itd*scat), 3E3*(ang+itd*scat), -50, 50 );
    trimadyvsty[itd] =
      TProfile( Form( "%smadyvsty", tds.c_str() ),
		Form( "%splet MAD(#Deltay) vs #theta_{y};%splet #theta_{y} [mrad];%splet MAD(#Deltay) [#mum]",
		      tds.c_str(), tds.c_str(), tds.c_str() ),
		100, -3E3*(ang+itd*scat), 3E3*(ang+itd*scat), 0, 50 );

    tridxCvsx[itd] = TProfile( Form( "%sdxCvsx", tds.c_str() ),
			       Form( "%splet dxC vs x;%splet x at C [mm];<%splets #Deltax C> [#mum]",
				     tds.c_str(), tds.c_str(), tds.c_str() ),
			       100, -midx[ipl], midx[ipl], -50, 50 );
    tridxCvsy[itd] = TProfile( Form( "%sdxCvsy", tds.c_str() ),
			       Form( "%splet dxC vs y;%splet y at C [mm];<%splets #Deltax C> [#mum]",
				     tds.c_str(), tds.c_str(), tds.c_str() ),
			       100, -midy[ipl], midy[ipl], -50, 50 );
    tridyCvsx[itd] = TProfile( Form( "%sdyCvsx", tds.c_str() ),
			       Form( "%splet dyC vs x;%splet x at C [mm];<%splets #Deltay C> [#mum]",
				     tds.c_str(), tds.c_str(), tds.c_str() ),
			       100, -midx[ipl], midx[ipl], -50, 50 );
    tridyCvsy[itd] = TProfile( Form( "%sdyCvsy", tds.c_str() ),
			       Form( "%splet dy vs y;%splet y at C [mm];<%splets #Deltay C> [#mum]",
				     tds.c_str(), tds.c_str(), tds.c_str() ),
			       100, -midy[ipl], midy[ipl], -50, 50 );
    tridxCvsax[itd] = TProfile( Form( "%sdxCvsax", tds.c_str() ),
				Form( "%splet dxC vs axAB;%splet AB slope x [mrad];<%splets #Deltax C> [#mum]",
				      tds.c_str(), tds.c_str(), tds.c_str() ),
				80, -2, 2, -50, 50 );
    tridyCvsay[itd] = TProfile( Form( "%sdyCvsay", tds.c_str() ),
				Form( "%splet dyC vs ayAB;%splet AB slope y [mrad];<%splets #Deltay C> [#mum]",
				      tds.c_str(), tds.c_str(), tds.c_str() ),
				80, -2, 2, -50, 50 );

    htrix[itd] = TH1I( Form( "%sx", tds.c_str() ),
		       Form( "%splet x;%splet x_{mid} [mm];%splets",
			     tds.c_str(), tds.c_str(), tds.c_str() ),
		       240, -12, 12 );
    htriy[itd] = TH1I( Form( "%sy", tds.c_str() ),
		       Form( "%splet y;%splet y_{mid} [mm];%splets",
			     tds.c_str(), tds.c_str(), tds.c_str() ),
		       120, -6, 6 );
    htrixy[itd] = new
      TH2I( Form( "%sxy", tds.c_str() ),
	    Form( "%splet x-y;%splet x_{mid} [mm];%splet y_{mid} [mm];%splets",
		  tds.c_str(), tds.c_str(), tds.c_str(), tds.c_str() ),
	    240, -12, 12, 120, -6, 6 );

    htritx[itd] = TH1I( Form( "%stx", tds.c_str() ),
			Form( "%splet #theta_{x};%splet #theta_{x} [mrad];%splets",
			      tds.c_str(), tds.c_str(), tds.c_str() ),
			200, -10E3*(ang+itd*scat), 10E3*(ang+itd*scat) );
    htrity[itd] = TH1I( Form( "%sty", tds.c_str() ),
			Form( "%splet #theta_{y};%splet #theta_{y} [mrad];%splets",
			      tds.c_str(), tds.c_str(), tds.c_str() ),
			200, -10E3*(ang+itd*scat), 10E3*(ang+itd*scat) );

    htrincol[itd] = TH1I( Form( "%sncol", tds.c_str() ),
			  Form( "%s cluster size x;columns/cluster;%s clusters on tracks",
				tds.c_str(), tds.c_str() ),
			  50, 0.5, 50.5 );
    htrinrow[itd] = TH1I( Form( "%snrow", tds.c_str() ),
			  Form( "%s cluster size y;rows/cluster;%s clusters on tracks",
				tds.c_str(), tds.c_str() ),
			  50, 0.5, 50.5 );
    htrinpix[itd] = TH1I( Form( "%snpix", tds.c_str() ),
			  Form( "%s cluster size;pixel/cluster;%s clusters on tracks",
				tds.c_str(), tds.c_str() ),
			  50, 0.5, 50.5 );

    hnpx1map[itd] = new
      TH2D( Form( "npx1map%i", ipl ),
	    Form( "plane %i 1-pix;x mod 18.4 [#mum];y mod 18.4 [#mum];plane %i 1-pix", ipl, ipl ),
	    74, 0, 1*ptchx[ipl]*1E3, 74, 0, 1*ptchy[ipl]*1E3 );
    hnpx2map[itd] = new
      TH2D( Form( "npx2map%i", ipl ),
	    Form( "plane %i 2-pix;x mod 18.4 [#mum];y mod 18.4 [#mum];plane %i 2-pix", ipl, ipl ),
	    74, 0, 1*ptchx[ipl]*1E3, 74, 0, 1*ptchy[ipl]*1E3 );
    hnpx3map[itd] = new
      TH2D( Form( "npx3map%i", ipl ),
	    Form( "plane %i 3-pix;x mod 18.4 [#mum];y mod 18.4 [#mum];plane %i 3-pix", ipl, ipl ),
	    74, 0, 1*ptchx[ipl]*1E3, 74, 0, 1*ptchy[ipl]*1E3 );
    hnpx4map[itd] = new
      TH2D( Form( "npx4map%i", ipl ),
	    Form( "plane %i 4-pix;x mod 18.4 [#mum];y mod 18.4 [#mum];plane %i 4-pix", ipl, ipl ),
	    74, 0, 1*ptchx[ipl]*1E3, 74, 0, 1*ptchy[ipl]*1E3 );

    trincolvsxm[itd] =
      TProfile( Form( "ncol%ivsxm", ipl ),
		Form( "plane %i cluster size;x mod 36.8 [#mum];plane %i <cluster size> [columns]", ipl, ipl ),
		74, 0, 2*ptchx[ipl]*1E3, 0, 2.5 );
    trinrowvsym[itd] =
      TProfile( Form( "nrow%ivsym", ipl ),
		Form( "plane %i cluster size;y mod 36.8 [#mum];plane %i <cluster size> [rows]", ipl, ipl ),
		74, 0, 2*ptchy[ipl]*1E3, 0, 2.5 );
    trinpixvsxmym[itd] = new
      TProfile2D( Form( "npix%ivsxmym", ipl ),
		  Form( "plane %i cluster size;x mod 73.6 [#mum];y mod 73.6 [#mum];plane %i <cluster size> [pixels]", ipl, ipl ),
		  74, 0, 4*ptchx[ipl]*1E3, 74, 0, 4*ptchy[ipl]*1E3, 0, 4.5 );
    trinpixgvsxmym[itd] = new
      TProfile2D( Form( "npix%igvsxmym", ipl ),
		  Form( "plane %i cluster size;x mod 73.6 [#mum];y mod 73.6 [#mum];plane %i <cluster size> [pixels]", ipl, ipl ),
		  74, 0, 4*ptchx[ipl]*1E3, 74, 0, 4*ptchy[ipl]*1E3, 0, 4.5 );

  } // triplets and driplets

  TH1I hntri( "ntri", "triplets per event;triplets;events", 51, -0.5, 50.5 );
  TH1I hndri( "ndri", "driplets per event;driplets;events", 51, -0.5, 50.5 );
  TProfile ntrivsev( "ntrivsev", "triplets per event vs time;time [events];<triplets/event>",
		     250, 0, 250*1000, -0.5, 99.5 );

  TProfile nclvspl( "nclvspl", "clusters per plane;plane;<clusters>", 6, 0.5, 6.5 );
  TProfile nlkvspl( "nlkvspl", "linked clusters per plane;plane;<linked clusters>", 6, 0.5, 6.5 );
  TProfile nfreevspl( "nfreevspl", "free clusters per plane;plane;<free clusters>", 6, 0.5, 6.5 );
  TH1I hnlkcl[6];
  for( unsigned jpl = 0; jpl < 6; ++jpl )
    hnlkcl[jpl] = TH1I( Form( "nlkcl%i", jpl+1 ),
			Form( "links per cluster in plane %i;links;plane %i clusters", jpl+1, jpl+1 ),
			10, -0.5,  9.5 );

  TH1I hexdx[9];
  TH1I hexdy[9];

  TH1I hexdxc[9];
  TH1I hexdyc[9];

  TProfile exdxvsy[9];
  TProfile exdyvsx[9];

  TProfile exdxvstx[9];
  TProfile exdyvsty[9];

  TProfile exmadxvstx[9];
  TProfile exmadyvsty[9];

  for( int ipl = 1; ipl <= 6; ++ipl ) {

    hexdx[ipl] = TH1I( Form( "exdx%i", ipl ),
		       Form( "ex dx %i;dx tri - plane %i [#mum];triplet - cluster pairs", ipl, ipl ),
		       200, -200*ff, 200*ff );
    hexdy[ipl] = TH1I( Form( "exdy%i", ipl ),
		       Form( "ex dy %i;dy tri - plane %i [#mum];triplet - cluster pairs", ipl, ipl ),
		       200, -200*ff, 200*ff );
    hexdxc[ipl] = TH1I( Form( "exdxc%i", ipl ),
			Form( "ex dx %i;dx tri - plane %i [#mum];triplet - cluster pairs", ipl, ipl ),
			200, -200*ff, 200*ff );
    hexdyc[ipl] = TH1I( Form( "exdyc%i", ipl ),
			Form( "ex dy %i;dy tri - plane %i [#mum];triplet - cluster pairs", ipl, ipl ),
			200, -200*ff, 200*ff );

    exdxvsy[ipl] = TProfile( Form( "exdxvsy%i", ipl ),
			     Form( "ex dx vs y %i;y at %i [mm];<#Deltax> [#mum]", ipl, ipl ),
			     100, -midy[ipl], midy[ipl], -200*ff, 200*ff );
    exdyvsx[ipl] = TProfile( Form( "exdyvsx%i", ipl ),
			     Form( "ex dy vs x %i;x at %i [mm];<#Deltay> [#mum]", ipl, ipl ),
			     100, -midx[ipl], midx[ipl], -200*ff, 200*ff );

    exdxvstx[ipl] =
      TProfile( Form( "exdxvstx%i", ipl ),
		Form( "dx vs tx at %i;slope x [mrad];<#Deltax> at %i [#mum]", ipl, ipl ),
		80, -3E3*ang, 3E3*ang, -100*ff, 100*ff );
    exdyvsty[ipl] =
      TProfile( Form( "exdyvsty%i", ipl ),
		Form( "dy vs ty at %i;slope y [mrad];<#Deltay> at %i [#mum]", ipl, ipl ),
		80, -3E3*ang, 3E3*ang, -100*ff, 100*ff );

    exmadxvstx[ipl] =
      TProfile( Form( "exmadxvstx%i", ipl ),
		Form( "MAD x vs tx at %i;slope x [mrad];MAD #Deltax at %i [#mum]", ipl, ipl ),
		80, -3E3*ang, 3E3*ang, 0, 100*ff );
    exmadyvsty[ipl] =
      TProfile( Form( "exmadyvsty%i", ipl ),
		Form( "MAD y vs ty at %i;slope y [mrad];MAD #Deltay at %i [#mum]", ipl, ipl ),
		80, -3E3*ang, 3E3*ang, 0, 100*ff );

  }  // ipl

  // driplet vertices:

  TH1I hdridd0( "dridd0",
		"driplets distance;driplets distance [mm];driplet pairs",
		200, -0.2, 0.2 );
  TH1I hdridd1( "dridd1",
		"driplets distance;driplets distance [mm];driplet pairs",
		200, -0.2, 0.2 );
  TH1I hdridd2( "dridd2",
		"driplets distance;driplets distance [mm];driplet pairs",
		200, -0.2, 0.2 );
  TH1I hdridd3( "dridd3",
		"driplets distance;driplets distance [mm];driplet pairs",
		200, -0.2, 0.2 );

  TH1I hdridz( "dridz",
	       "driplets intersection dz;driplets intersect dz [mm];driplet pairs",
	       200, -20, 20 );

  TH1I hdrizv( "drizv",
	       "driplets intersections;driplets intersect z [mm];driplet pairs",
	       5*int(zz[4]-zz[3]+20), zz[3], zz[4]+20 );

  TH1I hdridzm( "dridzm",
		"driplets mid intersection dz;driplets intersect dz [mm];mid driplet pairs",
		200, -20, 20 );

  TH1I htridxv( "tridxv",
		"triplet - vertex x;triplet - vertex x [mm];triplets",
		100, -0.5, 0.5 );
  TH1I htridyv( "tridyv",
		"triplet - vertex y;triplet - vertex y [mm];triplets",
		100, -0.5, 0.5 );

  TH1I hdridzt( "dridzt",
		"driplets mid intersection dz;driplets intersect dz [mm];mid driplet pairs",
		200, -20, 20 );

  TH1I hdrizvt( "drizvt",
		"driplets intersection with triplet;driplets intersect z [mm];driplet pairs with triplet",
		5*int(zz[4]-zz[3]+20), zz[3], zz[4]+20 );

  TH1I hdriddt( "driddt",
		"driplets distance;driplets distance [mm];driplet pairs",
		200, -0.2, 0.2 );
  TH1I hdrioat( "drioat",
		"driplet pair opening angle;opening angle [mrad];driplet pairs",
		200, 0, 100 );

  TH1I hdrizv0( "drizv0",
		"driplets intersection no triplet;driplets intersect z [mm];driplet pairs without triplet",
		5*int(zz[4]-zz[3]+20), zz[3], zz[4]+20 );

  TH1I hdridzd( "dridzd",
		"driplets mid intersection dz;driplets intersect dz [mm];mid driplet pairs",
		200, -20, 20 );

  TH1I hdrizvd( "drizvd",
		"driplets intersection;driplets intersect z [mm];driplet pairs",
		5*int(zz[4]-zz[3]+20), zz[3], zz[4]+20 );

  TH1I hdrioad( "drioad",
		"driplet pair opening angle;opening angle [mrad];driplet pairs",
		200, 0, 100 );

  TH1I hdrizvdl( "drizvdl",
		 "driplets intersection large angle;driplets intersect z [mm];driplet pairs",
		 5*int(zz[4]-zz[3]+20), zz[3], zz[4]+20 );
  TH1I hdrizvds( "drizvds",
		 "driplets intersection small angle;driplets intersect z [mm];driplet pairs",
		 5*int(zz[4]-zz[3]+20), zz[3], zz[4]+20 );

  TH1I hntrivert( "ntrivert", "triplets per driplet vertex;matched triplets;driplet vertices",
		  10, -0.5, 9.5 );

  TH1I hdrizvd0( "drizvd0",
		 "driplets intersection no triplet;driplets intersect z [mm];driplet pairs",
		 5*int(zz[4]-zz[3]+20), zz[3], zz[4]+20 );

  TH1I hdrizvdt( "drizvdt",
		 "driplets intersection with triplet;driplets intersect z [mm];driplet pairs",
		 5*int(zz[4]-zz[3]+20), zz[3], zz[4]+20 );

  TH1I hdridzdt( "dridzdt",
		 "driplets mid intersection dz;driplets intersect dz [mm];mid driplet pairs",
		 200, -20, 20 );
  TH1I hdrioadt( "drioadt",
		 "driplet pair opening angle;opening angle [mrad];driplet pairs",
		 200, 0, 100 );

  TH1I hdridxv( "dridxv",
		"driplet - vertex x;driplet - vertex x [mm];3rd driplets",
		100, -0.5, 0.5 );
  TH1I hdridyv( "dridyv",
		"driplet - vertex y;driplet - vertex y [mm];3rd driplets",
		100, -0.5, 0.5 );
  TH1I hndrivert( "ndrivert", "3rd driplets per vertex;matched extra driplets;driplet vertices",
		  10, -0.5, 9.5 );
  TH1I hdrizvdn( "drizvdn",
		 "driplets intersection no 3rd driplet;driplets intersect z [mm];driplet pairs",
		 5*int(zz[4]-zz[3]+20), zz[3], zz[4]+20 );
  TH1I hdrizvdd( "drizvdd",
		 "driplets intersection with extra driplet;driplets intersect z [mm];driplet pairs",
		 5*int(zz[4]-zz[3]+20), zz[3], zz[4]+20 );

  TH1I hdrizv01( "drizv01",
		 "driplets intersection;driplets intersect z [mm];driplet pairs",
		 5*int(zz[4]-zz[3]+20), zz[3], zz[4]+20 );
  TH1I hdrizv02( "drizv02",
		 "driplets intersection;driplets intersect z [mm];driplet pairs",
		 5*int(zz[4]-zz[3]+20), zz[3], zz[4]+20 );
  TH1I hdrizv05( "drizv05",
		 "driplets intersection;driplets intersect z [mm];driplet pairs",
		 5*int(zz[4]-zz[3]+20), zz[3], zz[4]+20 );
  TH1I hdrizv10( "drizv10",
		 "driplets intersection;driplets intersect z [mm];driplet pairs",
		 5*int(zz[4]-zz[3]+20), zz[3], zz[4]+20 );
  TH1I hdrizv20( "drizv20",
		 "driplets intersection;driplets intersect z [mm];driplet pairs",
		 5*int(zz[4]-zz[3]+20), zz[3], zz[4]+20 );
  TH1I hdrizv40( "drizv40",
		 "driplets intersection;driplets intersect z [mm];driplet pairs",
		 5*int(zz[4]-zz[3]+20), zz[3], zz[4]+20 );
  TH1I hdrizv80( "drizv80",
		 "driplets intersection;driplets intersect z [mm];driplet pairs",
		 5*int(zz[4]-zz[3]+20), zz[3], zz[4]+20 );

  TH1I hdrioaz( "drioaz", "driplet pair opening angle;opening angle [mrad];driplet pairs",
		200, 0, 100 );
  TH1I hvertx( "vertx", "driplet vertex x;vertex x [mm];driplet pairs", 240, -12, 12 );
  TH1I hverty( "verty", "driplet vertex y;vertex y [mm];driplet pairs", 160, -8, 8 );
  TH2I * hvertxy = new TH2I( "vertxy",
			     "driplet vertex map;vertex x [mm];vertex y [mm];driplet pairs",
			     120, -12, 12, 80, -8, 8 );

  TH1I hnvert( "nvert", "driplet vertices;vertices;events",
	       20, -0.5, 19.5 );

  // driplets - triplets

  TH1I hsixdx( "sixdx", "six dx;#Deltax [mm];triplet-driplet pairs", 800, -2*ff, 2*ff );
  TH1I hsixdy( "sixdy", "six dy;#Deltay [mm];triplet-driplet pairs", 400, -1*ff, 1*ff );

  TH1I hsixdxc( "sixdxc", "six dx;#Deltax [#mum];triplet-driplet pairs", 400, -5*sixcut*1E3, 5*sixcut*1E3 );
  TH1I hsixdyc( "sixdyc", "six dy;#Deltay [#mum];triplet-driplet pairs", 400, -5*sixcut*1E3, 5*sixcut*1E3 );

  TH2I * hsixxy = new
    TH2I( "sixxy", "sixplet x-y;sixplet x_{mid} [mm];sixplet y_{mid} [mm];sixplets",
	  240, -12, 12, 120, -6, 6 );

  TProfile sixdxvsx( "sixdxvsx",
		     "six #Deltax vs x;xB [mm];<driplet - triplet #Deltax [#mum]",
		     220, -11, 11, -100*ff, 100*ff );
  TProfile sixmadxvsx( "sixmadxvsx",
		       "six MAD x vs x;xB [mm];driplet - triplet MAD #Deltax [#mum]",
		       220, -11, 11, 0, 100*ff );
  TProfile sixmadxvsy( "sixmadxvsy",
		       "six MAD x vs y;yB [mm];driplet - triplet MAD #Deltax [#mum]",
		       106, -5.3, 5.3, 0, 100*ff );
  TProfile
    sixmadxvstx( "sixmadxvstx",
		 "six MAD x vs x;triplet #theta_{x} [mrad];driplet - triplet MAD #Deltax [#mum]",
		 80, -3E3*ang, 3E3*ang, 0, 100*ff );
  TProfile
    sixmadxvskx( "sixmadxvskx",
		 "six MAD x vs x;driplet-triplet #Delta#theta_{x} [mrad];driplet - triplet MAD #Deltax [#mum]",
		 80, -5E3*scat, 5E3*scat, 0, 100*ff );
  TProfile sixdxvsy( "sixdxvsy",
		     "six #Deltax vs y;yB [mm];<driplet - triplet #Deltax [#mum]",
		     106, -5.3, 5.3, -100*ff, 100*ff );
  TProfile sixdxvstx( "sixdxvstx",
		      "six #Deltax vs slope x;slope x [mrad];<driplet - triplet #Deltax> [#mum]",
		      80, -3E3*ang, 3E3*ang, -100*ff, 100*ff );

  TProfile sixdyvsx( "sixdyvsx",
		     "six #Deltay vs x;xB [mm];<driplet - triplet #Deltay [#mum]",
		     212, -10.6, 10.6, -100*ff, 100*ff );
  TProfile sixmadyvsx( "sixmadyvsx",
		       "six MAD y vs x;xB [mm];driplet - triplet MAD #Deltay [#mum]",
		       212, -10.6, 10.6, 0, 100*ff );

  TProfile sixdyvsy( "sixdyvsy",
		     "six #Deltay vs y;yB [mm];<driplet - triplet #Deltay [#mum]",
		     106, -5.3, 5.3, -100*ff, 100*ff );
  TProfile sixdyvsty( "sixdyvsty",
		      "six #Deltay vs slope y;slope y [mrad];<driplet - triplet #Deltay> [#mum]",
		      80, -3E3*ang, 3E3*ang, -100*ff, 100*ff );
  TProfile sixmadyvsy( "sixmadyvsy",
		       "six MAD y vs y;yB [mm];driplet - triplet MAD #Deltay [#mum]",
		       106, -5.3, 5.3, 0, 100*ff );
  TProfile
    sixmadyvsty( "sixmadyvsty",
		 "six MAD y vs #theta_{y};triplet #theta_{y} [mrad];driplet - triplet MAD #Deltay [#mum]",
		 80, -3E3*ang, 3E3*ang, 0, 100*ff );
  TProfile
    sixmadyvsky( "sixmadyvsky",
		 "six MAD y vs #Delta#theta_{y};driplet-triplet #Delta#theta_{y} [mrad];driplet - triplet MAD #Deltay [#mum]",
		 80, -2, 2, 0, 100*ff );

  TProfile2D * sixdxyvsxy = new
    TProfile2D( "sixdxyvsxy",
		"driplet - triplet #Delta_{xy} vs x-y;x_{mid} [mm];y_{mid} [mm];<sqrt(#Deltax^{2}+#Deltay^{2})> [mrad]",
		212, -10.6, 10.6, 106, -5.3, 5.3, 0, 100*ff );

  TH1I hsixkx( "sixkx",
	       "driplet slope x - triplet slope x;driplet slope x - triplet slope x [mrad];driplet-triplet pairs",
	       200, -20E3*scat, 20E3*scat );
  TH1I hsixky( "sixky",
	       "driplet slope y - triplet slope y;driplet slope y - triplet slope y [mrad];driplet-triplet pairs",
	       200, -20E3*scat, 20E3*scat );

  TProfile sixmadkxvsx( "sixmadkxvsx",
			"driplet - triplet slope_{x} vs x;x_{mid} [mm];MAD(#Delta#theta_{x}) [mrad]",
			424, -10.6, 10.6, 0, 20E3*scat );
  TProfile sixmadkyvsy( "sixmadkyvsy",
			"driplet - triplet slope_{y} vs y;y_{mid} [mm];MAD(#Delta#theta_{y}) [mrad]",
			212, -5.3, 5.3, 0, 20E3*scat );

  TProfile2D * sixdtvsxy = new
    TProfile2D( "sixdtvsxy",
		"driplet - triplet slope_{xy} vs x-y;x_{mid} [mm];y_{mid} [mm];<sqrt(#Delta#theta_{x}^{2}+#Delta#theta_{y}^{2})> [mrad]",
		424, -10.6, 10.6, 212, -5.3, 5.3, 0, 5E3*scat );

  TH1I hfourdx( "fourdx", "four dx;#Deltax [mm];doublet pairs", 800, -2*ff, 2*ff );
  TH1I hfourdy( "fourdy", "four dy;#Deltay [mm];doublet pairs", 400, -1*ff, 1*ff );

  TH1I hfourdxc( "fourdxc", "four dx;#Deltax [#mum];doublet pairs", 400, -5*sixcut*1E3, 5*sixcut*1E3 );
  TH1I hfourdyc( "fourdyc", "four dy;#Deltay [#mum];doublet pairs", 400, -5*sixcut*1E3, 5*sixcut*1E3 );

  TH1I hfourkx( "fourkx", "doublet kink x;doublet kink x [mrad];doublet pairs",
		500, -50E3*scat, 50E3*scat );
  TH1I hfourky( "fourky", "doublet kink y;doublet kink y [mrad];doublet pairs",
		500, -50E3*scat, 50E3*scat );
  TH2I * hfourkxky = new
    TH2I( "fourkxky", "doublet kink x-y;doublet kink x [mrad];doublet kink y [mrad];doublet pairs",
	  200, -20E3*scat, 20E3*scat, 200, -20E3*scat, 20E3*scat );

  TProfile fourkxvsx( "fourkxvsx",
		      "four kink_{x} vs x;x_{DUT} [mm];<#Delta#theta_{x}> [mrad]",
		      424, -10.6, 10.6 );
  TProfile fourkyvsy( "fourkyvsy",
		      "four kink_{y} vs y;y_{DUT} [mm];<#Delta#theta_{x}> [mrad]",
		      212, -5.3, 5.3 );
  TProfile fourmadkxvsx( "fourmadkxvsx",
			 "four MAD kink_{x} vs x;x_{DUT} [mm];MAD(#Delta#theta_{x}) [mrad]",
			 424, -10.6, 10.6 );
  TProfile fourmadkyvsx( "fourmadkyvsx",
			 "four MAD kink_{y} vs x;x_{DUT} [mm];MAD(#Delta#theta_{y}) [mrad]",
			 424, -10.6, 10.6 );
  TProfile fourmadkyvsy( "fourmadkyvsy",
			 "four MAD kink_{y} vs y;y_{DUT} [mm];MAD(#Delta#theta_{y}) [mrad]",
			 212, -5.3, 5.3 );

  TProfile2D * fourkvsxy = new
    TProfile2D( "fourkvsxy",
		"four kink_{xy} vs x-y;x_{DUT} [mm];y_{DUT} [mm];<sqrt(#Delta#theta_{x}^{2}+#Delta#theta_{y}^{2})> [mrad]",
		424, -10.6, 10.6, 212, -5.3, 5.3 );

  TH1I hfourdx0( "fourdx0", "four dx;#Deltax [#mum];doublet pairs", 400, -100, 100 );
  TH1I hfourdx5( "fourdx5", "four dx;#Deltax [#mum];doublet pairs", 400, -100, 100 );
  TH1I hfourdx3( "fourdx3", "four dx;#Deltax [#mum];doublet pairs", 400, -100, 100 );

  TH1I hfouroa( "fouroa", "four opening angle;opening angle [mrad];fours", 200, 0, 100 );
  TH1I hfouroa3( "fouroa3", "four opening angle at z3;opening angle [mrad];fours at z3", 200, 0, 100 );
  TH1I hfouroa4( "fouroa4", "four opening angle at z4;opening angle [mrad];fours at z4", 200, 0, 100 );

  TH2I * hfourzioa = new
    TH2I( "fourzioa",
	  "four intersection;four intersect z [mm];log_{10}(opening angle) [mrad];fours",
	  5*int(zz[4]-zz[3]+20), zz[3]-10, zz[4]+10, 30, -4, -1 );

  TH1I hfourzi0( "fourzi0", "four intersection;four intersect z [mm];fours",
		 5*int(zz[4]-zz[3]+20), zz[3]-10, zz[4]+10 );
  TH1I hfourzi1( "fourzi1", "four intersection;four intersect z [mm];fours",
		 5*int(zz[4]-zz[3]+20), zz[3]-10, zz[4]+10 );
  TH1I hfourzi2( "fourzi2", "four intersection;four intersect z [mm];fours",
		 5*int(zz[4]-zz[3]+20), zz[3]-10, zz[4]+10 );
  TH1I hfourzi4( "fourzi4", "four intersection;four intersect z [mm];fours",
		 5*int(zz[4]-zz[3]+20), zz[3]-10, zz[4]+10 );
  TH1I hfourzi6( "fourzi6", "four intersection;four intersect z [mm];fours",
		 5*int(zz[4]-zz[3]+20), zz[3]-10, zz[4]+10 );
  TH1I hfourzi8( "fourzi8", "four intersection;four intersect z [mm];fours",
		 5*int(zz[4]-zz[3]+20), zz[3]-10, zz[4]+10 );
  TH1I hfourzi12( "fourzi12", "four intersection;four intersect z [mm];fours",
		  5*int(zz[4]-zz[3]+20), zz[3]-10, zz[4]+10 );
  TH1I hfourzi16( "fourzi16", "four intersection;four intersect z [mm];fours",
		  5*int(zz[4]-zz[3]+20), zz[3]-10, zz[4]+10 );

  TH1I hndritri( "ndritri", "driplet matches per triplet;driplet matches;triplets",
		 9, -0.5, 8.5 );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // event loop:

  timespec ts;
  clock_gettime( CLOCK_REALTIME, &ts );
  time_t s0 = ts.tv_sec; // seconds since 1.1.1970
  long f0 = ts.tv_nsec; // nanoseconds
  double zeit1 = 0; // read
  double zeit2 = 0; // clus
  double zeit3 = 0; // track

  int iev = 0;
  uint64_t evTLU0 = 0;
  const double fTLU = 384E6; // 384 MHz TLU clock
  uint64_t prevTLU = 0;

  map < int, int > pxmap[9]; // for hot pixels

  list < vector < pixel > > pxlist[9];

  do {

    // Get next event:
    DetectorEvent evt = reader->GetDetectorEvent();

    if( evt.IsBORE() ) {
      cout << endl << "Begin Of Run Event" << endl << flush;
      eudaq::PluginManager::Initialize(evt);
    }

    if( lev < 99 )
      ldbg = 1;

    uint64_t evTLU = evt.GetTimestamp(); // 384 MHz = 2.6 ns
    if( iev < 2  ) // BORE has older time
      evTLU0 = evTLU;
    double evsec = (evTLU - evTLU0) / fTLU;
    t1Histo.Fill( evsec );
    t2Histo.Fill( evsec );
    t3Histo.Fill( evsec );
    t4Histo.Fill( evsec );
    t5Histo.Fill( evsec );

    double evdt = (evTLU - prevTLU) / fTLU;
    hdtus.Fill( evdt * 1E6 ); // [us]
    hdtms.Fill( evdt * 1E3 ); // [ms]
    prevTLU = evTLU;

    if( iev < 10 || ldbg )
      cout << "tele reading  " << run << "." << iev << "  taken " << evsec << endl;
    else if( iev < 100 && iev%10 == 0 )
      cout << "tele reading  " << run << "." << iev << "  taken " << evsec << endl;
    else if( iev < 1000 && iev%100 == 0 )
      cout << "tele reading  " << run << "." << iev << "  taken " << evsec << endl;
    else if( iev < 10000 && iev%1000 == 0 )
      cout << "tele reading  " << run << "." << iev << "  taken " << evsec << endl;
    else if( iev%10000 == 0 )
      cout << "tele reading  " << run << "." << iev << "  taken " << evsec << endl;

    StandardEvent sevt = eudaq::PluginManager::ConvertToStandard( evt );

    string MIM{"MIMOSA26"};
    int mpl = 1; // Mimosa planes start at 1

    if( ldbg ) cout << "planes " << sevt.NumPlanes() << endl << flush;

    for( size_t ipl = 0; ipl < sevt.NumPlanes(); ++ipl ) {

      const eudaq::StandardPlane &plane = sevt.GetPlane(ipl);

      if( ldbg )
	cout
	  << "  " << ipl
	  << ": plane " << plane.ID()
	  << " " << plane.Type() // NI
	  << " " << plane.Sensor() // MIMOSA26
	  << " frames " << plane.NumFrames() // 2
	  << " pivot " << plane.PivotPixel() // 6830
	  << " total " << plane.TotalPixels() // 663552
	  << " hits " << plane.HitPixels() // 5
	  ;

      //0: plane 1 NI MIMOSA26 frames 2 pivot 6830 total 663552 hits 5: 486 296 1: 635 68 1: 509 307 1: 510 307 1: 509 308 1

      if( plane.Sensor() != MIM ) {
	if( ldbg ) cout << endl;
	continue;
      }

      hpivot[mpl].Fill( plane.PivotPixel() );

      vector <double> pxl = plane.GetPixels<double>();

      hnpx0[mpl].Fill( pxl.size() );

      vector <pixel> pb; // for clustering

      for( size_t ipix = 0; ipix < pxl.size(); ++ipix ) {

	if( ldbg ) 
	  cout << " :"
	       << " " << plane.GetX(ipix) // col
	       << " " << plane.GetY(ipix) // row
	       << " " << plane.GetPixel(ipix) // charge
	       << flush;

	int ix = plane.GetX(ipix); // col pixel index
	int iy = plane.GetY(ipix); // row pixel index

	hcol0[mpl].Fill( ix );
	hrow0[mpl].Fill( iy );
	hbool[mpl].Fill( plane.GetPivot(ipix) );
	hmap0[mpl]->Fill( ix, iy );

	int ipx = ix * ny[mpl] + iy;

	if( ldbg )
	  cout << " " << ipx << flush;

	if( pxmap[mpl].count(ipx) )
	  ++pxmap[mpl][ipx];
	else
	  pxmap[mpl].insert( make_pair( ipx, 1 ) ); // slow

	if( hotset[mpl].count(ipx) ) {
	  if( ldbg )
	    cout << " hot" << flush;
	  continue; // skip hot
	}

	// fill pixel block for clustering:

	hcol[mpl].Fill( ix );
	hrow[mpl].Fill( iy );

	pixel px;
	px.col = ix;
	px.row = iy;
	px.nn = 1; // init neighbours for best resolution
	pb.push_back(px);

	if( pb.size() == 999 ) {
	  cout << "pixel buffer overflow in plane " << mpl
	       << ", event " << iev
	       << endl;
	  break;
	}

      } // pix

      if( ldbg )
	cout << " done" << endl << flush;
      
      hnpx[mpl].Fill( pb.size() );

      // for clustering:

      pxlist[mpl].push_back(pb);

      ++mpl;
      if( mpl > 6 ) break; // skip others

    } // planes

    ++iev;

  } while( reader->NextEvent() && iev < lev ); // event loop

  delete reader;

  clock_gettime( CLOCK_REALTIME, &ts );
  time_t s1 = ts.tv_sec; // seconds since 1.1.1970
  long f1 = ts.tv_nsec; // nanoseconds
  zeit1 += s1 - s0 + ( f1 - f0 ) * 1e-9; // read

  cout << "read " << iev << " events"
       << " in " << s1 - s0 + ( f1 - f0 ) * 1e-9 << " s"
       << endl;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // hot pixels:

  if( iev > 9999 ) {

    cout << endl << "Mimosa hot pixel list for run " << run << endl;

    ofstream hotFile( hotFileName.str() );

    hotFile << "# telescope hot pixel list for run " << run
	    << " with " << iev << " events"
	    << endl;

    for( int ipl = 1; ipl <= 6; ++ipl ) {
      hotFile << endl;
      hotFile << "plane " << ipl << endl;
      int nmax = 0;
      int ntot = 0;
      int nhot = 0;
      for( map < int, int >::iterator jpx = pxmap[ipl].begin(); jpx != pxmap[ipl].end(); ++ jpx ) {
	int nhit = jpx->second;
	ntot += nhit;
	if( nhit > nmax ) nmax = nhit;
	if( nhit > iev/128 ) {
	  ++nhot;
	  int ipx = jpx->first;
	  int ix = ipx/ny[ipl];
	  int iy = ipx%ny[ipl];
	  hotFile << "pix "
		  << setw(4) << ix
		  << setw(5) << iy
		  << "  " << nhit
		  << endl;
	}
      } // jpx
      cout
	<< "# " << ipl
	<< ": active " << pxmap[ipl].size()
	<< ", sum " << ntot
	<< ", max " << nmax
	<< ", hot " << nhot
	<< endl;
    } // ipl

    cout << "hot pixel list written to " << hotFileName.str() << endl;

    hotFile.close();

  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // make clusters:

  cout << endl << "parallel clustering" << flush;

  clock_gettime( CLOCK_REALTIME, &ts );
  time_t s2 = ts.tv_sec; // seconds since 1.1.1970
  long f2 = ts.tv_nsec; // nanoseconds

  list < vector <cluster> > clist[9];

  //#pragma omp sections // test, not parallel
#pragma omp parallel sections
  {
#pragma omp section
    {
      clist[1] = oneplane( 1, pxlist[1] );
    }
#pragma omp section
    {
      clist[2] = oneplane( 2, pxlist[2] );
    }
#pragma omp section
    {
      clist[3] = oneplane( 3, pxlist[3] );
    }
#pragma omp section
    {
      clist[4] = oneplane( 4, pxlist[4] );
    }
#pragma omp section
    {
      clist[5] = oneplane( 5, pxlist[5] );
    }
#pragma omp section
    {
      clist[6] = oneplane( 6, pxlist[6] );
    }

  } // parallel

  clock_gettime( CLOCK_REALTIME, &ts );
  time_t s3 = ts.tv_sec; // seconds since 1.1.1970
  long f3 = ts.tv_nsec; // nanoseconds
  zeit2 += s3 - s2 + ( f3 - f2 ) * 1e-9; // cluster

  cout << " in " << s3 - s2 + ( f3 - f2 ) * 1e-9 << " s" << endl;

  for( unsigned ipl = 0; ipl < 9; ++ipl )
    pxlist[ipl].clear(); // memory

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // alignment iterations:

  int maxiter = aligniteration + 1;
  if( maxiter < 9 ) maxiter = 9;

  for( ; aligniteration < maxiter; ++aligniteration ) {

    cout << endl << "using alignment iteration " << aligniteration << endl;

    for( unsigned ipl = 1; ipl <= 6; ++ipl ) {
      hxx[ipl]->Reset();
      hdx[ipl].Reset();
      hdy[ipl].Reset();
      dxvsx[ipl].Reset();
      dxvsy[ipl].Reset();
      dyvsx[ipl].Reset();
      dyvsy[ipl].Reset();
      hdx4[ipl].Reset();
      hdy4[ipl].Reset();
      dx4vsy[ipl].Reset();
      dy4vsx[ipl].Reset();
      effvsx[ipl].Reset();
      effvsy[ipl].Reset();
      effvsxm[ipl].Reset();
      effvsym[ipl].Reset();
      effvsxmym[ipl]->Reset();
    } // pl

    for( int itd = 0; itd < 2; ++itd ) { // triplets 0-1-2 and driplets 3-4-5

      hdxCA[itd].Reset();
      hdyCA[itd].Reset();
      dxCAvsx[itd].Reset();
      dyCAvsy[itd].Reset();

      htridx[itd].Reset();
      htridy[itd].Reset();

      htrikx[itd].Reset();
      htriky[itd].Reset();
      htrikxc[itd].Reset();
      htrikyc[itd].Reset();

      trimadkxvsx[itd].Reset();
      trimadkxvstx[itd].Reset();

      htridxc[itd].Reset();
      htridxci[itd].Reset();
      htridxct[itd].Reset();
      htridxctg[itd].Reset();

      tridxvsx[itd].Reset();
      tridxvsy[itd].Reset();
      tridxvsxm[itd].Reset();
      trimadxvsxm[itd].Reset();
      trimadxvsxmym[itd]->Reset();
      tridxvstx[itd].Reset();
      trimadxvstx[itd].Reset();

      htridxc1[itd].Reset();
      htridxc111[itd].Reset();
      htridxc2[itd].Reset();
      htridxc3[itd].Reset();
      htridxc4[itd].Reset();
      htridxc5[itd].Reset();
      htridxc6[itd].Reset();
      htridxc21[itd].Reset();
      htridxc22[itd].Reset();
      htridxc23[itd].Reset();
      htridxc24[itd].Reset();
      htridxc25[itd].Reset();
      htridxc223[itd].Reset();
      htridxc224[itd].Reset();
      htridxctgg[itd].Reset();

      htridyc[itd].Reset();
      htridyci[itd].Reset();
      htridyct[itd].Reset();
      htridyc1[itd].Reset();
      htridyc2[itd].Reset();
      htridyc3[itd].Reset();
      htridyc4[itd].Reset();
      htridyc5[itd].Reset();
      htridyc6[itd].Reset();
      htridycg[itd].Reset();

      tridyvsx[itd].Reset();
      tridyvsy[itd].Reset();
      tridyvsym[itd].Reset();
      trimadyvsym[itd].Reset();
      trimadyvsxmym[itd]->Reset();
      tridyvsty[itd].Reset();
      trimadyvsty[itd].Reset();

      htrix[itd].Reset();
      htriy[itd].Reset();
      htrixy[itd]->Reset();
      htritx[itd].Reset();
      htrity[itd].Reset();

      htrincol[itd].Reset();
      htrinrow[itd].Reset();
      htrinpix[itd].Reset();
      hnpx1map[itd]->Reset();
      hnpx2map[itd]->Reset();
      hnpx3map[itd]->Reset();
      hnpx4map[itd]->Reset();
      trincolvsxm[itd].Reset();
      trinrowvsym[itd].Reset();
      trinpixvsxmym[itd]->Reset();
      trinpixgvsxmym[itd]->Reset();

      tridxCvsx[itd].Reset();
      tridxCvsy[itd].Reset();
      tridyCvsx[itd].Reset();
      tridyCvsy[itd].Reset();
      tridxCvsax[itd].Reset();
      tridyCvsay[itd].Reset();

    } // itd

    hntri.Reset();
    ntrivsev.Reset();
    hndri.Reset();

    for( unsigned ipl = 1; ipl <= 6; ++ipl ) {

      hexdx[ipl].Reset();
      hexdy[ipl].Reset();
      hexdxc[ipl].Reset();
      exdxvsy[ipl].Reset();
      exdxvstx[ipl].Reset();
      exmadxvstx[ipl].Reset();
      hexdyc[ipl].Reset();
      exdyvsx[ipl].Reset();
      exdyvsty[ipl].Reset();
      exmadyvsty[ipl].Reset();

    }

    hdridd0.Reset();
    hdridd1.Reset();
    hdridd2.Reset();
    hdridd3.Reset();

    hdridz.Reset();
    hdrizv.Reset();
    hdridzm.Reset();

    htridxv.Reset();
    htridyv.Reset();

    hdridzt.Reset();
    hdrizvt.Reset();
    hdriddt.Reset();
    hdrioat.Reset();
    hdrizv0.Reset();

    hdridzd.Reset();
    hdrizvd.Reset();
    hdrioad.Reset();
    hdrizvdl.Reset();
    hdrizvds.Reset();

    hntrivert.Reset();

    hdrizvd0.Reset();

    hdridzdt.Reset();
    hdrizvdt.Reset();
    hdrioadt.Reset();

    hdridxv.Reset();
    hdridyv.Reset();
    hndrivert.Reset();

    hdrizvdn.Reset();
    hdrizvdd.Reset();

    hdrizv01.Reset();
    hdrizv02.Reset();
    hdrizv05.Reset();
    hdrizv10.Reset();
    hdrizv20.Reset();
    hdrizv40.Reset();
    hdrizv80.Reset();

    hdrioaz.Reset();
    hvertx.Reset();
    hverty.Reset();
    hvertxy->Reset();

    hnvert.Reset();

    hsixdx.Reset();
    hsixdy.Reset();
    hsixdxc.Reset();
    sixdxvsx.Reset();
    sixmadxvsx.Reset();
    sixdxvsy.Reset();
    sixdxvstx.Reset();
    sixmadxvsy.Reset();
    sixmadxvstx.Reset();
    sixmadxvskx.Reset();

    hsixdyc.Reset();
    sixdyvsx.Reset();
    sixmadyvsx.Reset();
    sixdyvsy.Reset();
    sixdyvsty.Reset();
    sixmadyvsy.Reset();
    sixmadyvsty.Reset();
    sixmadyvsky.Reset();

    hsixxy->Reset();
    sixdxyvsxy->Reset();

    sixmadkxvsx.Reset();
    sixmadkyvsy.Reset();
    hsixkx.Reset();
    hsixky.Reset();
    sixdtvsxy->Reset();

    hfourdx.Reset();
    hfourdy.Reset();
    hfourdxc.Reset();
    hfourdyc.Reset();
    hfourkx.Reset();
    hfourkxky->Reset();
    fourkxvsx.Reset();
    fourkyvsy.Reset();
    fourkvsxy->Reset();
    fourmadkxvsx.Reset();
    hfourdx0.Reset();
    hfourdx5.Reset();
    hfourdx3.Reset();
    hfourky.Reset();
    fourmadkyvsx.Reset();

    hfourzi0.Reset();
    hfourzi1.Reset();
    hfourzi2.Reset();
    hfourzi4.Reset();
    hfourzi6.Reset();
    hfourzi8.Reset();
    hfourzi12.Reset();
    hfourzi16.Reset();
    hfouroa.Reset();
    hfouroa3.Reset();
    hfouroa4.Reset();
    hfourzioa->Reset();

    hndritri.Reset();

    // loop over events, correlate planes:

    list < vector <cluster> >::iterator evi[9];
    for( unsigned ipl = 1; ipl <= 6; ++ipl )
      evi[ipl] = clist[ipl].begin();

    iev = 0;

    cout << "tracking ev";

    for( ; evi[1] != clist[1].end();
	 ++evi[1], ++evi[2], ++evi[3], ++evi[4], ++evi[5], ++evi[6] ) {

      vector <cluster> cl[9]; // Mimosa planes
      for( unsigned ipl = 1; ipl <= 6; ++ipl )
	cl[ipl] = *evi[ipl];

      ++iev;
      if( iev%10000 == 0 )
	cout << " " << iev << flush;

      // final cluster plots:

      if( aligniteration == maxiter-1 ) {

	for( unsigned ipl = 1; ipl <= 6; ++ipl ) {

	  hncl[ipl].Fill( cl[ipl].size() );

	  for( vector<cluster>::iterator cA = cl[ipl].begin(); cA != cl[ipl].end(); ++cA ) {

	    hccol[ipl].Fill( cA->col );
	    hcrow[ipl].Fill( cA->row );
	    unsigned nrow = cA->scr/(1024*1024);
	    unsigned ncol = (cA->scr - nrow*1024*1024)/1024;
	    unsigned npix = cA->scr % 1024;
	    hnpix[ipl].Fill( npix );
	    hncol[ipl].Fill( ncol );
	    hnrow[ipl].Fill( nrow );
	    hmindxy[ipl].Fill( cA->mindxy );

	  } // clus

	} // planes

      } // laster iter

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // cluster pair correlations:

      double twocut = 2; // [mm]
      if( aligniteration >= 1 )
	twocut = 0.5;

      for( int itd = 0; itd < 2; ++itd ) { // triplets 0-1-2 and driplets 3-4-5

	int im = 2; // mid plane triplet
	int ibeg = 1;
	int iend = 3;
	if( itd == 1 ) {
	  im = 5; // mid plane driplet
	  ibeg = 4;
	  iend = 6;
	}

	// A = mid plane:

	for( vector<cluster>::iterator cA = cl[im].begin(); cA != cl[im].end(); ++cA ) {

	  double xA = cA->col*ptchx[im] - alignx[im];
	  double yA = cA->row*ptchy[im] - aligny[im];
	  double xmid = xA - midx[im];
	  double ymid = yA - midy[im];
	  xA = xmid - ymid*rotx[im];
	  yA = ymid + xmid*roty[im];

	  for( int ipl = ibeg; ipl <= iend; ++ipl ) {

	    if( ipl == im ) continue;

	    double sign = ipl - im; // along track: -1 or 1

	    // B = A +- 1

	    for( vector<cluster>::iterator cB = cl[ipl].begin(); cB != cl[ipl].end(); ++cB ) {

	      double xB = cB->col*ptchx[ipl] - alignx[ipl];
	      double yB = cB->row*ptchy[ipl] - aligny[ipl];
	      double xmid = xB - midx[ipl];
	      double ymid = yB - midy[ipl];
	      xB = xmid - ymid*rotx[ipl];
	      yB = ymid + xmid*roty[ipl];

	      double dx = xB - xA;
	      double dy = yB - yA;
	      hxx[ipl]->Fill( xA, xB );
	      if( fabs(dy) < twocut ) {
		hdx[ipl].Fill( dx ); // for shift: fixed sign
		dxvsx[ipl].Fill( xB, dx*sign ); // for turn angle: sign along track
		dxvsy[ipl].Fill( yB, dx      );
	      }
	      if( fabs(dx) < twocut ) {
		hdy[ipl].Fill( dy );
		dyvsx[ipl].Fill( xB, dy      );
		dyvsy[ipl].Fill( yB, dy*sign ); // for tilt angle: sign along track
	      }

	    } // clusters

	  } // ipl

	} // cl mid

      } // upstream and downstream internal correlations

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // telescope plane efficiency:

      double isoCut = 0.30/ptchx[1]; // [px]
      double effCut = 0.100*ff; // [mm]

      for( int ipl = 1; ipl <= 6; ++ipl ) {

	int ib = 2;
	int im = 3; // mid plane triplet
	int ie = 4;
	if(      ipl == 2 ) {
	  ib = 1;
	  im = 3;
	  ie = 4;
	}
	else if( ipl == 3 ) {
	  ib = 1;
	  im = 2;
	  ie = 4;
	}
	else if( ipl == 4 ) {
	  ib = 2;
	  im = 3;
	  ie = 5;
	}
	else if( ipl == 5 ) {
	  ib = 3;
	  im = 4;
	  ie = 6;
	}
	else if( ipl == 6 ) {
	  ib = 3;
	  im = 4;
	  ie = 5;
	}

	double zD = zz[ipl] + alignz[ipl];

	for( vector<cluster>::iterator cA = cl[ib].begin(); cA != cl[ib].end(); ++cA ) {

	  if( cA->mindxy < isoCut ) continue;

	  double xA = cA->col*ptchx[ib] - alignx[ib];
	  double yA = cA->row*ptchy[ib] - aligny[ib];
	  double zA = zz[ib] + alignz[ib];
	  double xmid = xA - midx[ib];
	  double ymid = yA - midy[ib];
	  xA = xmid - ymid*rotx[ib];
	  yA = ymid + xmid*roty[ib];

	  for( vector<cluster>::iterator cC = cl[ie].begin(); cC != cl[ie].end(); ++cC ) {

	    if( cC->mindxy < isoCut ) continue;

	    double xC = cC->col*ptchx[ie] - alignx[ie];
	    double yC = cC->row*ptchy[ie] - aligny[ie];
	    double zC = zz[ie] + alignz[ie];
	    double xmid = xC - midx[ie];
	    double ymid = yC - midy[ie];
	    xC = xmid - ymid*rotx[ie];
	    yC = ymid + xmid*roty[ie];

	    double dx2 = xC - xA;
	    double dy2 = yC - yA;
	    double dzCA = zC - zA;

	    if( fabs( dx2 ) > ang * dzCA ) continue; // angle cut
	    if( fabs( dy2 ) > ang * dzCA ) continue; // angle cut

	    double xavg2 = 0.5*(xA + xC);
	    double yavg2 = 0.5*(yA + yC);
	    double zavg2 = 0.5*(zA + zC);

	    double slpx = ( xC - xA ) / dzCA; // slope x
	    double slpy = ( yC - yA ) / dzCA; // slope y

	    for( vector<cluster>::iterator cB = cl[im].begin(); cB != cl[im].end(); ++cB ) {

	      double xB = cB->col*ptchx[im] - alignx[im]; // stretch and shift
	      double yB = cB->row*ptchy[im] - aligny[im];
	      double zB = zz[im] + alignz[im];
	      double xmid = xB - midx[im];
	      double ymid = yB - midy[im];
	      xB = xmid - ymid*rotx[im];
	      yB = ymid + xmid*roty[im];

	      // kinks:

	      double kx = (xC-xB)/(zC-zB) - (xB-xA)/(zB-zA);
	      double ky = (yC-yB)/(zC-zB) - (yB-yA)/(zB-zA);

	      if( fabs( kx ) > tricut ) continue;
	      if( fabs( ky ) > tricut ) continue;

	      // inter/extrapolate track to D:

	      double da = zD - zavg2;
	      double xi = xavg2 + slpx * da; // triplet at D
	      double yi = yavg2 + slpy * da;

	      // transform into local frame:

	      double xr = xi + yi*rotx[ipl] + alignx[ipl];
	      double yr = yi - xi*roty[ipl] + aligny[ipl];

	      if( fabs( xr ) > 10.4 ) continue; // fiducial
	      if( fabs( yr ) >  5.2 ) continue; // fiducial

	      // eff pl:

	      int nm = 0;

	      for( vector<cluster>::iterator cD = cl[ipl].begin(); cD != cl[ipl].end(); ++cD ) {

		double xD = cD->col*ptchx[ipl] - midx[ipl];
		double yD = cD->row*ptchy[ipl] - midy[ipl];

		double dx4 = xD - xr;
		double dy4 = yD - yr;

		if( fabs( dy4 ) < effCut ) {
		  hdx4[ipl].Fill( dx4*1E3 );
		  dx4vsy[ipl].Fill( yr, dx4*1E3 );
		}
		if( fabs( dx4 ) < effCut ) {
		  hdy4[ipl].Fill( dy4*1E3 );
		  dy4vsx[ipl].Fill( xr, dy4*1E3 );
		}

		if( fabs( dx4 ) > effCut ) continue;
		if( fabs( dy4 ) > effCut ) continue;

		++nm;

		if( nm > 0 ) break; // one link is enough

	      } // cl D

	      effvsx[ipl].Fill( xr, nm );
	      effvsy[ipl].Fill( yr, nm );

	      double xmod2 = fmod( xr + sizex[ipl] + 0.5*ptchx[ipl], 2*ptchx[ipl] );
	      double ymod2 = fmod( yr + sizey[ipl] + 0.5*ptchy[ipl], 2*ptchy[ipl] );
	      effvsxm[ipl].Fill( xmod2*1E3, nm );
	      effvsym[ipl].Fill( ymod2*1E3, nm );
	      effvsxmym[ipl]->Fill( xmod2*1E3, ymod2*1E3, nm );

	    } // cl B

	  } // cl C

	} // cl A

      } // eff planes

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // triplets 2 vs 3-1:
      // driplets 5 vs 6-4:

      vector <triplet> triplets;
      vector <triplet> driplets;

      for( int itd = 0; itd < 2; ++itd ) { // triplets 0-1-2 and driplets 3-4-5

	int ib = 1;
	int ie = 3;
	int im = 2; // mid plane triplet

	if( itd == 1 ) {
	  ib = 4;
	  ie = 6;
	  im = 5; // mid plane driplet
	}

	double zA = zz[ib] + alignz[ib];
	double zC = zz[ie] + alignz[ie];
	double zB = zz[im] + alignz[im];
	double dzCA = zC - zA;

	for( vector<cluster>::iterator cA = cl[ib].begin(); cA != cl[ib].end(); ++cA ) {

	  double xA = cA->col*ptchx[ib] - alignx[ib];
	  double yA = cA->row*ptchy[ib] - aligny[ib];
	  double xmid = xA - midx[ib];
	  double ymid = yA - midy[ib];
	  xA = xmid - ymid*rotx[ib];
	  yA = ymid + xmid*roty[ib];

	  unsigned nrowA = cA->scr/(1024*1024);
	  unsigned ncolA = (cA->scr - nrowA*1024*1024)/1024;
	  //unsigned npixA = cA->scr % 1024;
	  bool goodncolA = 1;
	  if( ncolA > 4 ) goodncolA = 0;
	  if( ncolA == 2 && nrowA < 3 ) goodncolA = 0;
	  bool goodnrowA = 1;
	  if( nrowA > 4 ) goodnrowA = 0;
	  if( nrowA == 2 && nrowA < 3 ) goodnrowA = 0;

	  for( vector<cluster>::iterator cC = cl[ie].begin(); cC != cl[ie].end(); ++cC ) {

	    double xC = cC->col*ptchx[ie] - alignx[ie];
	    double yC = cC->row*ptchy[ie] - aligny[ie];
	    double xmid = xC - midx[ie];
	    double ymid = yC - midy[ie];
	    xC = xmid - ymid*rotx[ie];
	    yC = ymid + xmid*roty[ie];

	    double dx2 = xC - xA;
	    double dy2 = yC - yA;

	    hdxCA[itd].Fill( dx2 );
	    hdyCA[itd].Fill( dy2 );

	    if( fabs( dy2 ) < 0.001 * dzCA )
	      dxCAvsx[itd].Fill( xC, dx2 );

	    if( fabs( dx2 ) < 0.001 * dzCA )
	      dyCAvsy[itd].Fill( yC, dy2 );

	    if( fabs( dx2 ) > ( 5*ang + 50*itd*scat ) * dzCA ) continue; // beam divergence and scat
	    if( fabs( dy2 ) > ( 5*ang + 50*itd*scat ) * dzCA ) continue;

	    double xavg2 = 0.5*(xA + xC);
	    double yavg2 = 0.5*(yA + yC);
	    double zavg2 = 0.5*(zA + zC);

	    double slpx = ( xC - xA ) / dzCA; // slope x
	    double slpy = ( yC - yA ) / dzCA; // slope y

	    // interpolate track to B:

	    double dz = zB - zavg2;
	    double xm = xavg2 + slpx * dz; // triplet at B
	    double ym = yavg2 + slpy * dz;

	    // transform into local frame:

	    double xr = xm + ym*rotx[im] + alignx[im];
	    double yr = ym - xm*roty[im] + aligny[im];

	    double xmod1 = fmod( xr + sizex[im] + 0.5*ptchx[im], 1*ptchx[im] );
	    double ymod1 = fmod( yr + sizey[im] + 0.5*ptchy[im], 1*ptchy[im] );
	    double xmod2 = fmod( xr + sizex[im] + 0.5*ptchx[im], 2*ptchx[im] );
	    double ymod2 = fmod( yr + sizey[im] + 0.5*ptchy[im], 2*ptchy[im] );
	    double xmod4 = fmod( xr + sizex[im] + 0.5*ptchx[im], 4*ptchx[im] );
	    double ymod4 = fmod( yr + sizey[im] + 0.5*ptchy[im], 4*ptchy[im] );

	    unsigned nrowC = cC->scr/(1024*1024);
	    unsigned ncolC = (cC->scr - nrowC*1024*1024)/1024;
	    //unsigned npixC = cC->scr % 1024;
	    bool goodncolC = 1;
	    if( ncolC > 4 ) goodncolC = 0;
	    if( ncolC == 2 && nrowC < 3 ) goodncolC = 0;
	    bool goodnrowC = 1;
	    if( nrowC > 4 ) goodnrowC = 0;
	    if( nrowC == 2 && nrowC < 3 ) goodnrowC = 0;

	    for( vector<cluster>::iterator cB = cl[im].begin(); cB != cl[im].end(); ++cB ) {

	      double xB = cB->col*ptchx[im] - alignx[im]; // stretch and shift
	      double yB = cB->row*ptchy[im] - aligny[im];
	      double xmid = xB - midx[im];
	      double ymid = yB - midy[im];
	      xB = xmid - ymid*rotx[im];
	      yB = ymid + xmid*roty[im];

	      double dxm = xB - xm;
	      double dym = yB - ym;

	      htridx[itd].Fill( dxm*1E3 );
	      htridy[itd].Fill( dym*1E3 );

	      // kinks:

	      double kx = (xC-xB)/(zC-zB) - (xB-xA)/(zB-zA);
	      double ky = (yC-yB)/(zC-zB) - (yB-yA)/(zB-zA);

	      htrikx[itd].Fill( kx*1E3 );
	      htriky[itd].Fill( ky*1E3 );

	      bool iso = 1;
	      if( cA->mindxy < isoCut ) iso = 0;
	      if( cC->mindxy < isoCut ) iso = 0;
	      if( cB->mindxy < isoCut ) iso = 0;

	      unsigned nrowB = cB->scr/(1024*1024);
	      unsigned ncolB = (cB->scr - nrowB*1024*1024)/1024;
	      unsigned npixB = cB->scr % 1024;

	      if( ncolB > 99 )
		cout << "scrB " << cB->scr
		     << ", nrow " << nrowB
		     << ", ncol " << ncolB
		     << ", npix " << npixB
		     << endl;

	      bool goodncolB = 1;
	      if( ncolB > 4 ) goodncolB = 0;
	      if( ncolB == 2 && nrowB < 3 ) goodncolB = 0;
	      bool goodnrowB = 1;
	      if( nrowB > 4 ) goodnrowB = 0;
	      if( nrowB == 2 && nrowB < 3 ) goodnrowB = 0;

	      if( fabs( ky ) < tricut ) {

		htrikxc[itd].Fill( kx*1E3 );
		if( iso )
		  htrikxci[itd].Fill( kx*1E3 );

		trimadkxvsx[itd].Fill( xr, fabs(kx)*1E3 ); // flat
		trimadkxvstx[itd].Fill( slpx*1E3, fabs(kx)*1E3 ); // U-shape

		htridxc[itd].Fill( dxm*1E3 );
		if( iso )
		  htridxci[itd].Fill( dxm*1E3 );

		tridxvsx[itd].Fill( xr, dxm*1E3 );
		tridxvsy[itd].Fill( yr, dxm*1E3 );
		tridxvstx[itd].Fill( slpx*1E3, dxm*1E3 ); // adjust zpos, same sign

		tridxvsxm[itd].Fill( xmod2*1E3, dxm*1E3 );

		trimadxvsxm[itd].Fill( xmod2*1E3, fabs(dxm)*1E3 );
		trimadxvsxmym[itd]->Fill( xmod2*1E3, ymod2*1E3, fabs(dxm)*1E3 );
		trimadxvstx[itd].Fill( slpx*1E3, fabs(dxm)*1E3 ); // U-shape

		if( fabs( slpx ) < ang ) {

		  htridxct[itd].Fill( dxm*1E3 ); // 3.95

		  if( goodncolA && goodncolC ) { // position bias?

		    htridxctg[itd].Fill( dxm*1E3 ); // 3.70 (24%)

		    if( goodncolB )
		      htridxctgg[itd].Fill( dxm*1E3 ); // 3.13 thr 4, 3.70 thr 5

		  } // good A, C

		  htrixm[itd].Fill( xmod1*1E3 );

		  if(      ncolB == 1 ) {
		    htridxc1[itd].Fill( dxm*1E3 ); // 2.40
		    htrixm1[itd].Fill( xmod1*1E3 );
		  }
		  else if( ncolB == 2 ) {
		    htridxc2[itd].Fill( dxm*1E3 ); // 4.14
		    htrixm2[itd].Fill( xmod1*1E3 );
		  }
		  else if( ncolB == 3 ) {
		    htridxc3[itd].Fill( dxm*1E3 ); // 3.30
		    htrixm3[itd].Fill( xmod1*1E3 );
		  }
		  else if( ncolB == 4 ) {
		    htridxc4[itd].Fill( dxm*1E3 ); // 3.23
		    htrixm4[itd].Fill( xmod1*1E3 );
		  }
		  else if( ncolB == 5 )
		    htridxc5[itd].Fill( dxm*1E3 ); // 5.59
		  else
		    htridxc6[itd].Fill( dxm*1E3 ); // 

		  if( ncolB == 2 ) {

		    if(      nrowB == 1 )
		      htridxc21[itd].Fill( dxm*1E3 ); // 4.93
		    else if( nrowB == 2 )
		      htridxc22[itd].Fill( dxm*1E3 ); // 4.26 most
		    else if( nrowB == 3 )
		      htridxc23[itd].Fill( dxm*1E3 ); // 3.27
		    else if( nrowB == 4 )
		      htridxc24[itd].Fill( dxm*1E3 ); // 2.35
		    else
		      htridxc25[itd].Fill( dxm*1E3 ); // 

		    if( nrowB == 2 ) {
		      if( npixB == 3 )
			htridxc223[itd].Fill( dxm*1E3 ); // 4.16
		      else
			htridxc224[itd].Fill( dxm*1E3 ); // 4.32 most
		    }

		  }

		  if( ncolA == 1 && ncolB == 1 &&  ncolC == 1 )
		    htridxc111[itd].Fill( dxm*1E3 ); // 1.7

		} // slpx

	      } // dy

	      if( fabs( kx ) < tricut ) {

		htrikyc[itd].Fill( ky*1E3 );
		if( iso )
		  htrikyci[itd].Fill( ky*1E3 );
		htridyc[itd].Fill( dym*1E3 );
		if( iso )
		  htridyci[itd].Fill( dym*1E3 );
		tridyvsx[itd].Fill( xr, dym*1E3 );
		tridyvsy[itd].Fill( yr, dym*1E3 );
		tridyvsym[itd].Fill( ymod2*1E3, dym*1E3 );
		trimadyvsym[itd].Fill( ymod2*1E3, fabs(dym)*1E3 );
		trimadyvsxmym[itd]->Fill( xmod2*1E3, ymod2*1E3, fabs(dym)*1E3 );
		tridyvsty[itd].Fill( slpy*1E3, dym*1E3 );
		trimadyvsty[itd].Fill( slpy*1E3, fabs(dym)*1E3 ); // U-shape
		if( fabs( slpy ) < 0.001 )
		  htridyct[itd].Fill( dym*1E3 );

		if( fabs( slpy ) < 0.001 ) {

		  if(      nrowB == 1 )
		    htridyc1[itd].Fill( dym*1E3 ); // 
		  else if( nrowB == 2 )
		    htridyc2[itd].Fill( dym*1E3 ); // 
		  else if( nrowB == 3 )
		    htridyc3[itd].Fill( dym*1E3 ); // 
		  else if( nrowB == 4 )
		    htridyc4[itd].Fill( dym*1E3 ); // 
		  else if( nrowB == 5 )
		    htridyc5[itd].Fill( dym*1E3 ); // 
		  else
		    htridyc6[itd].Fill( dym*1E3 ); // 

		  if( goodnrowA && goodnrowB && goodnrowC )
		    htridycg[itd].Fill( dym*1E3 ); // 

		} // slpy

	      } // dx

	      // cut x and y:

	      if( fabs( kx ) < tricut &&
		  fabs( ky ) < tricut ) {

		htrincol[itd].Fill( ncolB );
		htrinrow[itd].Fill( nrowB );
		htrinpix[itd].Fill( npixB );

		if(      npixB == 1 )
		  hnpx1map[itd]->Fill( xmod1*1E3, ymod1*1E3 );
		else if( npixB == 2 )
		  hnpx2map[itd]->Fill( xmod1*1E3, ymod1*1E3 );
		else if( npixB == 3 )
		  hnpx3map[itd]->Fill( xmod1*1E3, ymod1*1E3 );
		else if( npixB == 4 )
		  hnpx4map[itd]->Fill( xmod1*1E3, ymod1*1E3 );

		if( fabs( slpx ) < 0.001 )
		  trincolvsxm[itd].Fill( xmod2*1E3, ncolB );

		if( fabs( slpy ) < 0.001 )
		  trinrowvsym[itd].Fill( ymod2*1E3, nrowB );

		if( fabs( slpx ) < 0.001 && fabs( slpy ) < 0.001 )
		  trinpixvsxmym[itd]->Fill( xmod4*1E3, ymod4*1E3, npixB );

		if( goodnrowA && goodnrowC && goodncolA && goodncolC ) // better resolution
		  trinpixgvsxmym[itd]->Fill( xmod4*1E3, ymod4*1E3, npixB );

		// store triplets:

		triplet tri;
		tri.xm = xavg2;
		tri.ym = yavg2;
		tri.zm = zavg2;
		tri.sx = slpx;
		tri.sy = slpy;
		tri.rxy = sqrt(dxm*dxm+dym*dym);
		tri.kx = kx;
		tri.ky = ky;
		tri.ghost = 0;

		tri.i1 = distance( cl[ib].begin(), cA ); // starts a 0
		tri.i2 = distance( cl[im].begin(), cB );
		tri.i3 = distance( cl[ie].begin(), cC );

		vector <double> ux(3);
		ux[0] = xA;
		ux[1] = xB;
		ux[2] = xC;
		tri.vx = ux;

		vector <double> uy(3);
		uy[0] = yA;
		uy[1] = yB;
		uy[2] = yC;
		tri.vy = uy;

		if( itd )
		  driplets.push_back(tri);
		else
		  triplets.push_back(tri);

		htrix[itd].Fill( xavg2 );
		htriy[itd].Fill( yavg2 );
		htrixy[itd]->Fill( xavg2, yavg2 );
		htritx[itd].Fill( slpx*1E3 );
		htrity[itd].Fill( slpy*1E3 );

		// check z spacing: A-B as baseline

		double dzAB = zB - zA;
		double ax = ( xB - xA ) / dzAB; // slope x
		double ay = ( yB - yA ) / dzAB; // slope y
		double dz = zC - zB;
		double xk = xB + ax * dz; // at C
		double yk = yB + ay * dz; // at C
		double dx = xC - xk;
		double dy = yC - yk;
		tridxCvsx[itd].Fill( xk, dx*1E3 );
		tridxCvsy[itd].Fill( yk, dx*1E3 );
		tridyCvsx[itd].Fill( xk, dy*1E3 );
		tridyCvsy[itd].Fill( yk, dy*1E3 );
		tridxCvsax[itd].Fill( ax*1E3, dx*1E3 ); // adjust zpos, same sign
		tridyCvsay[itd].Fill( ay*1E3, dy*1E3 );

	      } // triplet

	    } // cl B

	  } // cl C

	} // cl A

      } // triplets and driplets

      // eliminate ghosts:

      for( unsigned jB = 0; jB < driplets.size(); ++jB ) { // j = B = downstream

	for( unsigned jj = jB + 1; jj < driplets.size(); ++jj ) {

	  //for( unsigned ipl = 0; ipl < 3; ++ipl ) {
	  for( unsigned ipl = 2; ipl < 3; ++ipl ) {

	    if( pow( driplets[jj].vx[ipl] - driplets[jB].vx[ipl], 2 ) +
		pow( driplets[jj].vy[ipl] - driplets[jB].vy[ipl], 2 ) < 4*ptchx[4]*ptchy[4] ) {
	      if( driplets[jj].rxy > driplets[jB].rxy )
		driplets[jj].ghost = 1;
	      else
		driplets[jB].ghost = 1;
	    }

	  } // ipl

	} // jj

      } // jB

      hntri.Fill( triplets.size() );
      ntrivsev.Fill( iev, triplets.size() );
      hndri.Fill( driplets.size() );

      if( lev < 100 )
	cout << " " << iev
	     << " cl " << cl[1].size()
	     << "  " << cl[2].size()
	     << "  " << cl[3].size()
	     << "  " << cl[4].size()
	     << "  " << cl[5].size()
	     << "  " << cl[6].size()
	     << ", tri " << triplets.size()
	     << ", dri " << driplets.size()
	     << endl;

      // cluster usage in tracks:

      for( unsigned ipl = 1; ipl <= 6; ++ipl )
	nclvspl.Fill( ipl, cl[ipl].size() );

      vector < vector <int> > nlkpl(6); // per plane, starting at 0
      for( unsigned jpl = 0; jpl < 6; ++jpl ) {
	nlkpl[jpl] = vector <int> ( cl[jpl+1].size() );
      }

      for( unsigned int iA = 0; iA < triplets.size(); ++iA ) { // 1,2,3
	++nlkpl[0].at( triplets[iA].i1 ); // tracks per cluster
	++nlkpl[1].at( triplets[iA].i2 );
	++nlkpl[2].at( triplets[iA].i3 );
      }

      for( unsigned int jB = 0; jB < driplets.size(); ++jB ) { // 4,5,6
	++nlkpl[3].at( driplets[jB].i1 );
	++nlkpl[4].at( driplets[jB].i2 );
	++nlkpl[5].at( driplets[jB].i3 );
      }

      for( unsigned jpl = 0; jpl < 6; ++jpl ) {

	int nlk = 0;
	int nfree = 0;

	for( unsigned icl = 0; icl < cl[jpl+1].size(); ++icl ) {

	  if( nlkpl[jpl].at(icl) )
	    ++nlk;
	  else
	    ++nfree;

	  hnlkcl[jpl+1].Fill( nlkpl[jpl].at(icl) );

	} // icl

	nlkvspl.Fill( jpl+1, nlk );
	nfreevspl.Fill( jpl+1, nfree );

      } // jpl from 0

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // extrapolate triplets to each downstream plane
      // dy vs ty: dz

      for( unsigned int iA = 0; iA < triplets.size(); ++iA ) { // i = A = upstream

	double xmA = triplets[iA].xm;
	double ymA = triplets[iA].ym;
	double zmA = triplets[iA].zm;
	double sxA = triplets[iA].sx;
	double syA = triplets[iA].sy;

	for( int ipl = 4; ipl <= 6; ++ipl ) {

	  // triplet at plane:

	  double zA = zz[ipl] + alignz[ipl] - zmA; // z from mid of triplet to plane
	  double xA = xmA + sxA * zA; // triplet at mid
	  double yA = ymA + syA * zA;

	  for( vector<cluster>::iterator cC = cl[ipl].begin(); cC != cl[ipl].end(); ++cC ) {

	    double xC = cC->col*ptchx[ipl] - alignx[ipl];
	    double yC = cC->row*ptchy[ipl] - aligny[ipl];
	    double xmid = xC - midx[ipl];
	    double ymid = yC - midy[ipl];
	    xC = xmid - ymid*rotx[ipl];
	    yC = ymid + xmid*roty[ipl];

	    double dx = xC - xA;
	    double dy = yC - yA;
	    hexdx[ipl].Fill( dx*1E3 );
	    hexdy[ipl].Fill( dy*1E3 );
	    if( fabs( dy ) < 0.5 ) {
	      hexdxc[ipl].Fill( dx*1E3 );
	      exdxvsy[ipl].Fill( yC, dx*1E3 );
	      exdxvstx[ipl].Fill( sxA*1E3, dx*1E3 );
	      exmadxvstx[ipl].Fill( sxA*1E3, fabs(dx)*1E3 );
	    }
	    if( fabs( dx ) < 0.5 ) {
	      hexdyc[ipl].Fill( dy*1E3 );
	      exdyvsx[ipl].Fill( xC, dy*1E3 );
	      exdyvsty[ipl].Fill( syA*1E3, dy*1E3 );
	      exmadyvsty[ipl].Fill( syA*1E3, fabs(dy)*1E3 );
	    }

	  } // clus

	} // planes

      } // triplets

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // extrapolate driplets to each upstream plane

      vector <vertex> vv;

      for( unsigned int jB = 0; jB < driplets.size(); ++jB ) { // j = B = downstream

	if( driplets[jB].ghost ) continue;

	double xmB = driplets[jB].xm;
	double ymB = driplets[jB].ym;
	double zmB = driplets[jB].zm;
	double sxB = driplets[jB].sx;
	double syB = driplets[jB].sy;

	for( int ipl = 1; ipl <= 3; ++ipl ) {

	  // driplet at plane:

	  double zB = zz[ipl] + alignz[ipl] - zmB; // z from mid of driplet to plane
	  double xB = xmB + sxB * zB; // driplet at mid
	  double yB = ymB + syB * zB;

	  for( vector<cluster>::iterator cC = cl[ipl].begin(); cC != cl[ipl].end(); ++cC ) {

	    double xC = cC->col*ptchx[ipl] - alignx[ipl];
	    double yC = cC->row*ptchy[ipl] - aligny[ipl];
	    double xmid = xC - midx[ipl];
	    double ymid = yC - midy[ipl];
	    xC = xmid - ymid*rotx[ipl];
	    yC = ymid + xmid*roty[ipl];

	    double dx = xC - xB;
	    double dy = yC - yB;
	    hexdx[ipl].Fill( dx*1E3 );
	    hexdy[ipl].Fill( dy*1E3 );
	    if( fabs( dy ) < 0.5 ) {
	      hexdxc[ipl].Fill( dx*1E3 );
	      exdxvsy[ipl].Fill( yC, dx*1E3 );
	      exdxvstx[ipl].Fill( sxB*1E3, dx*1E3 );
	      exmadxvstx[ipl].Fill( sxB*1E3, fabs(dx)*1E3 );
	    }
	    if( fabs( dx ) < 0.5 ) {
	      hexdyc[ipl].Fill( dy*1E3 );
	      exdyvsx[ipl].Fill( xC, dy*1E3 );
	      exdyvsty[ipl].Fill( syB*1E3, dy*1E3 );
	      exmadyvsty[ipl].Fill( syB*1E3, fabs(dy)*1E3 );
	    }

	  } // clus

	} // planes

	// vertices:

	for( unsigned int jj = jB + 1; jj < driplets.size(); ++jj ) {

	  if( driplets[jj].ghost ) continue;

	  if( pow( driplets[jj].vx[0] - driplets[jB].vx[0], 2 ) +
	      pow( driplets[jj].vy[0] - driplets[jB].vy[0], 2 ) < 4*ptchx[4]*ptchy[4] ) continue;

	  if( pow( driplets[jj].vx[1] - driplets[jB].vx[1], 2 ) +
	      pow( driplets[jj].vy[1] - driplets[jB].vy[1], 2 ) < 4*ptchx[4]*ptchy[4] ) continue;

	  // opening angles:

	  double sxj = driplets[jj].sx;
	  double syj = driplets[jj].sy;

	  if( fabs( sxj - sxB ) < 0.002 ) continue;
	  if( fabs( syj - syB ) < 0.002 ) continue;

	  double xmj = driplets[jj].xm;
	  double ymj = driplets[jj].ym;
	  double zmj = driplets[jj].zm;

	  // g = (xm,ym,zm) + dz*(sx,sy,1)

	  // cross product:

	  double Nx = syj*1 - 1*syB;
	  double Ny = 1*sxB - sxj*1;
	  double Nz = sxj*syB - syj*sxB;
	  double Nn = sqrt( Nx*Nx + Ny*Ny + Nz*Nz );

	  // distance:

	  double dd = ( (xmj-xmB)*Nx + (ymj-ymB)*Ny ) / Nn;

	  hdridd0.Fill( dd );

	  hdridd1.Fill( dd );

	  // z intersection:

	  double oax = sxj - sxB;
	  double oay = syj - syB;
	  double zix = ( xmB - xmj ) / oax;
	  double ziy = ( ymB - ymj ) / oay;

	  double dzi = zix - ziy;

	  hdridz.Fill( dzi ); // broad

	  // weight by opening angle:

	  double zv = zmB + ( oax*oax*zix + oay*oay*ziy ) / (oax*oax+oay*oay);

	  hdrizv.Fill( zv );

	  if( zv > zz[3]+alignz[3] && zv < zz[4]+alignz[4] ) {
	    hdridzm.Fill( dzi );
	    hdridd2.Fill( dd );
	  }

	  // opening angle 3D:

	  double snB = sqrt( sxB*sxB + syB*syB + 1 );
	  double snj = sqrt( sxj*sxj + syj*syj + 1 );
	  double sprod = sxB*sxj + syB*syj + 1;
	  double oa = acos( sprod / snB/snj );

	  // vertex xy:

	  double lzB = zv - zmB; // z from mid of driplet to intersect
	  double xvB = xmB + sxB * lzB; // vertex
	  double yvB = ymB + syB * lzB;

	  double lzj = zv - zmj; // z from mid of driplet to intersect
	  double xvj = xmj + sxj * lzj; // vertex
	  double yvj = ymj + syj * lzj;

	  double xv = 0.5 * ( xvB + xvj );
	  double yv = 0.5 * ( yvB + yvj );
	    
	  // incoming triplet:

	  int nmt = 0;

	  for( unsigned int iA = 0; iA < triplets.size(); ++iA ) { // i = A = upstream

	    double xmA = triplets[iA].xm;
	    double ymA = triplets[iA].ym;
	    double zmA = triplets[iA].zm;
	    double sxA = triplets[iA].sx;
	    double syA = triplets[iA].sy;

	    // triplet at vertex:

	    double lz = zv - zmA; // z from mid of triplet to vertex
	    double xA = xmA + sxA * lz; // triplet at vertex
	    double yA = ymA + syA * lz;

	    double dxv = xA - xv;
	    double dyv = yA - yv;

	    if( fabs(dyv) < 0.1 )
	      htridxv.Fill( dxv );
	    if( fabs(dxv) < 0.1 )
	      htridyv.Fill( dyv );

	    if( fabs(dxv) < 0.060 && fabs(dyv) < 0.060 )
	      ++nmt;

	  } // iA

	  if( nmt ) {
	    hdridzt.Fill( dzi );
	    hdrizvt.Fill( zv );
	    hdriddt.Fill( dd );
	    hdrioat.Fill( oa*1E3 ); // [mrad]
	  }
	  else
	    hdrizv0.Fill( zv ); // ghosts?

	  if( fabs(dd) < 0.050 ) {

	    hdridzd.Fill( dzi );
	    hdrizvd.Fill( zv ); // material peaks
	    hdrioad.Fill( oa*1E3 ); // [mrad]

	    if( oa > 0.017 ) // 1 deg
	      hdrizvdl.Fill( zv ); // material peaks
	    else
	      hdrizvds.Fill( zv ); // material peaks

	    hntrivert.Fill( nmt );

	    if( nmt ) {
	      hdrizvdt.Fill( zv ); // material peaks
	      hdridzdt.Fill( dzi );
	      hdrioadt.Fill( oa*1E3 ); // [mrad] mean 18 at 2.4 GeV with 3 mm Pb
	      hvertx.Fill( xv );
	      hverty.Fill( yv );
	      hvertxy->Fill( xv, yv );
	    }
	    else
	      hdrizvd0.Fill( zv ); // ghosts?

	    vertex vx;
	    vx.x = xv;
	    vx.y = yv;
	    vx.z = zv;
	    vx.i = jB;
	    vx.j = jj;
	    vv.push_back( vx );
	    
	    // 3rd driplet = scattered e

	    int nmd = 0;

	    for( unsigned int kk = 0; kk < driplets.size(); ++kk ) {

	      if( driplets[kk].ghost ) continue;

	      if( kk == jB ) continue;
	      if( kk == jj ) continue;

	      double xmk = driplets[kk].xm;
	      double ymk = driplets[kk].ym;
	      double zmk = driplets[kk].zm;
	      double sxk = driplets[kk].sx;
	      double syk = driplets[kk].sy;

	      // driplet at vertex:

	      double lz = zv - zmk; // z from mid of driplet to vertex
	      double xk = xmk + sxk * lz; // driplet at vertex
	      double yk = ymk + syk * lz;

	      double dxv = xk - xv;
	      double dyv = yk - yv;

	      hdridxv.Fill( dxv );
	      hdridyv.Fill( dyv );

	      if( fabs(dxv) < 0.060 && fabs(dyv) < 0.060 )
		++nmd;

	    } // iA

	    hndrivert.Fill( nmd );

	    if( nmd ) {
	      hdrizvdd.Fill( zv );
	    }
	    else
	      hdrizvdn.Fill( zv );

	  } // cut dd

	  if( fabs(dzi) < 0.1 )
	    hdrizv01.Fill( zv ); // material peaks

	  if( fabs(dzi) < 0.2 )
	    hdrizv02.Fill( zv ); // material peaks

	  if( fabs(dzi) < 0.5 )
	    hdrizv05.Fill( zv ); // material peaks

	  if( fabs(dzi) < 1.0 )
	    hdrizv10.Fill( zv ); // material peaks

	  if( fabs(dzi) < 2.0 )
	    hdrizv20.Fill( zv ); // material peaks

	  if( fabs(dzi) < 4.0 )
	    hdrizv40.Fill( zv ); // material peaks

	  if( fabs(dzi) < 8.0 )
	    hdrizv80.Fill( zv ); // material peaks

	  if( zv > zz[3]+alignz[3] && zv < zz[4]+alignz[4] && fabs(dzi) < 4 ) {

	    hdridd3.Fill( dd ); // clean

	    hdrioaz.Fill( oa*1E3 ); // [mrad]

	  } // z-vertex

	} // jj

      } // driplets

      hnvert.Fill( vv.size() );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // match triplets and driplets, measure kink angles
      // beware of parallaxe: match triplet and driplet at the z intersect

      for( unsigned int iA = 0; iA < triplets.size(); ++iA ) { // i = A = upstream

	double xmA = triplets[iA].xm;
	double ymA = triplets[iA].ym;
	double zmA = triplets[iA].zm;
	double sxA = triplets[iA].sx;
	double syA = triplets[iA].sy;

	// re-define triplet to doublet: less telescope material

	double sxt = ( triplets[iA].vx[2] - triplets[iA].vx[1] ) / dz32;
	double syt = ( triplets[iA].vy[2] - triplets[iA].vy[1] ) / dz32;
	double xmt = ( triplets[iA].vx[2] + triplets[iA].vx[1] ) * 0.5;
	double ymt = ( triplets[iA].vy[2] + triplets[iA].vy[1] ) * 0.5;

	vector <int> vjB;

	for( unsigned int jB = 0; jB < driplets.size(); ++jB ) { // j = B = downstream

	  if( driplets[jB].ghost ) continue;

	  double xmB = driplets[jB].xm;
	  double ymB = driplets[jB].ym;
	  double zmB = driplets[jB].zm;
	  double sxB = driplets[jB].sx;
	  double syB = driplets[jB].sy;

	  // 2-plane driplet:

	  double sxd = ( driplets[jB].vx[1] - driplets[jB].vx[0] ) / dz54;
	  double syd = ( driplets[jB].vy[1] - driplets[jB].vy[0] ) / dz54;
	  double xmd = ( driplets[jB].vx[1] + driplets[jB].vx[0] ) * 0.5;
	  double ymd = ( driplets[jB].vy[1] + driplets[jB].vy[0] ) * 0.5;

	  // opening angles:

	  double oax = sxd - sxt;
	  double oay = syd - syt;

	  double snt = sqrt( sxt*sxt + syt*syt + 1 );
	  double snd = sqrt( sxd*sxd + syd*syd + 1 );
	  double sprod = sxt*sxd + syt*syd + 1;
	  double oa = acos( sprod / snt/snd ); // positive

	  // triplet at mid driplet:

	  double wz = zm45 - zm23;
	  double xe = xmt + sxt * wz;
	  double ye = ymt + syt * wz;

	  // z intersection:

	  double zix = DUTz-zm45;
	  double ziy = DUTz-zm45;
	  if( fabs(oax) > 1e-4 )
	    zix = ( xe - xmd ) / oax;
	  if( fabs(oay) > 1e-4 )
	    ziy = ( ye - ymd ) / oay;

	  // weight by opening angle:

	  double zi = zm45 + ( oax*oax*zix + oay*oay*ziy ) / (oax*oax+oay*oay);
	  if( zi > zm45 ) zi = DUTz;
	  if( zi < zm23 ) zi = DUTz;

	  // driplet at intersect:

	  double lz = zi - zmB; // z from mid of triplet to mid
	  double xB = xmB + sxB * lz; // triplet at mid
	  double yB = ymB + syB * lz;

	  lz = zi - zm45;
	  double xd = xmd + sxd * lz;
	  double yd = ymd + syd * lz;

	  // triplet at DUT:

	  lz = zi - zmA; // z from mid of triplet to mid driplet
	  double xA = xmA + sxA * lz; // triplet at mid
	  double yA = ymA + syA * lz;

	  lz = zi - zm23;
	  double xt = xmt + sxt * lz;
	  double yt = ymt + syt * lz;

	  double dx4 = xd-xt;
	  double dy4 = yd-yt;

	  // driplet - triplet:

	  double dx = xB - xA;
	  double dy = yB - yA;
	  double dxy = sqrt( dx*dx + dy*dy );

	  double kx = sxB - sxA;
	  double ky = syB - syA;
	  double kxy = sqrt( kx*kx + ky*ky );

	  hsixdx.Fill( dx ); // [mm] for align fit
	  hsixdy.Fill( dy ); // [mm] for align fit

	  if( fabs(dy) < sixcut ) {

	    hsixdxc.Fill( dx*1E3 );

	    sixdxvsx.Fill( xA, dx*1E3 );
	    sixmadxvsx.Fill( xA, fabs(dx)*1E3 );
	    sixdxvsy.Fill( yA, dx*1E3 );
	    sixdxvstx.Fill( sxA*1E3, dx*1E3 );
	    sixmadxvsy.Fill( yA, fabs(dx)*1E3 );
	    sixmadxvstx.Fill( sxA*1E3, fabs(dx)*1E3 );
	    sixmadxvskx.Fill( kx*1E3, fabs(dx)*1E3 ); // U-shape

	    if( yA > -4.6 && yA < 5.1 )
	      sixmadkxvsx.Fill( xA, fabs(kx)*1E3 );

	    if( xA > -8.6 && xA < 7.4 &&
		yA > -4.6 && yA < 5.1 ) {
	      hsixkx.Fill( kx*1E3 );
	      hsixky.Fill( ky*1E3 );
	    }

	  } // dy

	  if( fabs(dx) < sixcut ) {

	    hsixdyc.Fill( dy*1E3 );

	    sixdyvsx.Fill( xA, dy*1E3 );
	    sixmadyvsx.Fill( xA, fabs(dy)*1E3 );
	    sixdyvsy.Fill( yA, dy*1E3 );
	    sixdyvsty.Fill( syA*1E3, dy*1E3 );
	    sixmadyvsy.Fill( yA, fabs(dy)*1E3 );
	    sixmadyvsty.Fill( syA*1E3, fabs(dy)*1E3 );
	    sixmadyvsky.Fill( ky*1E3, fabs(dy)*1E3 ); // U-shape

	    if( xA > -8.6 && xA < 7.4 )
	      sixmadkyvsy.Fill( yA, fabs(ky)*1E3 );

	    if( xA > -8.6 && xA < 7.4 &&
		yA > -4.6 && yA < 5.1 ) {
	      hsixkx.Fill( kx*1E3 );
	      hsixky.Fill( ky*1E3 );
	    }

	  }

	  // match:

	  if( fabs(dy) < sixcut &&
	      fabs(dx) < sixcut ) {

	    hsixxy->Fill( xA, yA );
	    sixdxyvsxy->Fill( xA, yA, dxy );

	    sixdtvsxy->Fill( xA, yA, kxy*1E3 );

	    vjB.push_back(jB); // brems conversions: looser sixcut?

	  } // match

	  hfourdx.Fill( dx4 );
	  hfourdy.Fill( dy4 );

	  double kx4 = sxd - sxt;
	  double ky4 = syd - syt;
	  double kxy4 = sqrt( kx4*kx4 + ky4*ky4 );

	  if( fabs(dy4) < sixcut ) {

	    hfourdxc.Fill( dx4*1E3 );

	    if( fabs(triplets[iA].kx) < 3*scatm && fabs(driplets[jB].kx) < 3*scatm )
	      hfourdx0.Fill( dx4*1E3 );
	    else if( fabs(triplets[iA].kx) > 5*scatm || fabs(driplets[jB].kx) > 5*scatm )
	      hfourdx5.Fill( dx4*1E3 );
	    else
	      hfourdx3.Fill( dx4*1E3 );

	  }

	  if( fabs(dx4) < sixcut ) {
	    hfourdyc.Fill( dy4*1E3 );
	  }

	  if( fabs(dx4) < sixcut &&fabs(dy4) < sixcut ) {

	    fourkxvsx.Fill( xt, kx4*1E3 );
	    fourkyvsy.Fill( yt, ky4*1E3 );
	    fourkvsxy->Fill( xt, yt, kxy4*1E3 );
	    fourmadkxvsx.Fill( xt, fabs(kx4)*1E3 );
	    fourmadkyvsy.Fill( yt, fabs(ky4)*1E3 );

	    if( xA > -10.6 + 3*trng && xA < 10.6 - 3*trng )
	      hfourkx.Fill( kx4*1E3 );

	    if( yA > -5.3 + 3*trng && yA < 5.3 - 3*trng ) {
	      hfourky.Fill( ky4*1E3 );
	      fourmadkyvsx.Fill( xt, fabs(ky4)*1E3 );
	    }
	    if( xA > -10.6 + 3*trng && xA < 10.6 - 3*trng &&
		yA > -5.3 + 3*trng && yA < 5.3 - 3*trng )
	      hfourkxky->Fill( kx4*1E3, ky4*1E3 );

	    if( fabs(oax) > 1e-4 && fabs(oay) > 1e-4 ) {

	      hfourzi0.Fill( zi );

	      hfouroa.Fill( oa*1E3 );

	      if( zi > z3-2 && zi < z3+5 )
		hfouroa3.Fill( oa*1E3 );
	      if( zi > z4-2 && zi < z4+2 )
		hfouroa4.Fill( oa*1E3 );

	      hfourzioa->Fill( zi, log(oa)/log10 );

	      if( oa > 0.001 )
		hfourzi1.Fill( zi );
	      if( oa > 0.002 )
		hfourzi2.Fill( zi );
	      if( oa > 0.004 )
		hfourzi4.Fill( zi );
	      if( oa > 0.006 )
		hfourzi6.Fill( zi );
	      if( oa > 0.008 )
		hfourzi8.Fill( zi );
	      if( oa > 0.012 )
		hfourzi12.Fill( zi );
	      if( oa > 0.016 )
		hfourzi16.Fill( zi );

	    } // oax, oay

	  } // match

	} // driplets

	hndritri.Fill( vjB.size() );

      } // triplets

    } // events

    cout << endl;

    clock_gettime( CLOCK_REALTIME, &ts );
    time_t s3 = ts.tv_sec; // seconds since 1.1.1970
    long f3 = ts.tv_nsec; // nanoseconds
    zeit3 += s3 - s2 + ( f3 - f2 ) * 1e-9; // track

    clock_gettime( CLOCK_REALTIME, &ts );
    time_t s9 = ts.tv_sec; // seconds since 1.1.1970
    long f9 = ts.tv_nsec; // nanoseconds

    cout << endl
	 << "done after " << iev << " events"
	 << " in " << s9 - s0 + ( f9 - f0 ) * 1e-9 << " s"
	 << " (read " << zeit1
	 << " s, cluster " << zeit2
	 << " s, tracking " << zeit3 << " s)"
	 << endl;

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // alignment fits:

    if( maxiter == aligniteration + 1 ) {      
      histoFile->Write(); // before fitting
      histoFile->Close();
    }

    cout << endl << "alignment fits:" << endl;

    if( aligniteration < 3 ) {

      for( int ipl = 1; ipl <= 6; ++ipl ) {

	double nb = hdx[ipl].GetNbinsX();
	double ne = hdx[ipl].GetSumOfWeights();
	double nm = hdx[ipl].GetMaximum();

	cout << endl << hdx[ipl].GetTitle() << endl;
	if( nm < 99 ) {
	  cout << "  peak " << nm << " not enough" << endl;
	  continue;
	}

	double fwhmx = GetFWHM( hdx[ipl] );
	double xm = hdx[ipl].GetBinCenter( hdx[ipl].GetMaximumBin() );
	cout << "  Inside  " << ne << " (" << ne/nb << " per bin)" << endl;
	cout << "  Maximum " << nm << " (factor " << nm/ne*nb << " above mean)" << endl;
	cout << "  at " << xm << endl;
	cout << "  fwhm " << fwhmx << endl;

	TF1 * fgp0x = new TF1( "fgp0x", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", xm-fwhmx, xm+fwhmx );
	fgp0x->SetParameter( 0, nm ); // amplitude
	fgp0x->SetParameter( 1, xm );
	fgp0x->SetParameter( 2, 0.5*fwhmx ); // sigma
	fgp0x->SetParameter( 3, hdx[ipl].GetBinContent(1) ); // BG
	hdx[ipl].Fit( "fgp0x", "qr" );
	cout << "Fit Gauss + BG:"
	     << endl << "  A " << fgp0x->GetParameter(0)
	     << endl << "mid " << fgp0x->GetParameter(1) << " mm"
	     << endl << "sig " << fgp0x->GetParameter(2) << " mm"
	     << endl << " BG " << fgp0x->GetParameter(3)
	     << endl;

	alignx[ipl] += fgp0x->GetParameter(1);

	delete fgp0x;

	// dy:

	nb = hdy[ipl].GetNbinsX();
	ne = hdy[ipl].GetSumOfWeights();
	nm = hdy[ipl].GetMaximum();
	cout << endl << hdy[ipl].GetTitle() << endl;
	double fwhmy = GetFWHM( hdy[ipl] );
	double ym = hdy[ipl].GetBinCenter( hdy[ipl].GetMaximumBin() );
	cout << "  Inside  " << ne << " (" << ne/nb << " per bin)" << endl;
	cout << "  Maximum " << nm << " (factor " << nm/ne*nb << " above mean)" << endl;
	cout << "  at " << ym << endl;
	cout << "  fwhm " << fwhmy << endl;

	TF1 * fgp0y = new TF1( "fgp0y", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", ym-fwhmy, ym+fwhmy );
	fgp0y->SetParameter( 0, nm ); // amplitude
	fgp0y->SetParameter( 1, ym );
	fgp0y->SetParameter( 2, 0.5*fwhmy ); // sigma
	fgp0y->SetParameter( 3, hdy[ipl].GetBinContent(1) ); // BG
	hdy[ipl].Fit( "fgp0y", "qr" );
	cout << "Fit Gauss + BG:"
	     << endl << "  A " << fgp0y->GetParameter(0)
	     << endl << "mid " << fgp0y->GetParameter(1) << " mm"
	     << endl << "sig " << fgp0y->GetParameter(2) << " mm"
	     << endl << " BG " << fgp0y->GetParameter(3)
	     << endl;

	aligny[ipl] += fgp0y->GetParameter(1);

	delete fgp0y;

	// x-y rotation:

	if( aligniteration >= 1 && dxvsy[ipl].GetEntries() > 999 ) {

	  dxvsy[ipl].Fit( "pol1", "q", "", -midy[ipl], midy[ipl] );
	  TF1 * fdxvsy = dxvsy[ipl].GetFunction( "pol1" );
	  cout << endl << dxvsy[ipl].GetTitle()
	       << " slope " << fdxvsy->GetParameter(1)
	       << " = rotx"
	       << endl;
	  rotx[ipl] += fdxvsy->GetParameter(1);

	}

	if( aligniteration >= 1 && dyvsx[ipl].GetEntries() > 999 ) {

	  dyvsx[ipl].Fit( "pol1", "q", "", -midx[ipl], midx[ipl] );
	  TF1 * fdyvsx = dyvsx[ipl].GetFunction( "pol1" );
	  cout << endl << dyvsx[ipl].GetTitle()
	       << " slope " << fdyvsx->GetParameter(1)
	       << " = roty"
	       << endl;
	  roty[ipl] -= fdyvsx->GetParameter(1); // sign

	}

      } // ipl

    } // iter

    // tridx:

    if( aligniteration >= 3 && htridx[1].GetMaximum() > 99 ) {

      for( int itd = 0; itd < 2; ++itd ) { // triplets 0-1-2 and driplets 3-4-5

	int ib = 1;
	int im = 2; // mid plane triplet
	int ie = 3;
	if( itd == 1 ) {
	  ib = 4;
	  im = 5; // mid plane driplet
	  ie = 6;
	}

	double nb = htridxc[itd].GetNbinsX();
	double ne = htridxc[itd].GetSumOfWeights();
	double nm = htridxc[itd].GetMaximum();
	double fwhmx = GetFWHM( htridxc[itd] );
	double xm = htridxc[itd].GetBinCenter( htridxc[itd].GetMaximumBin() );
	cout << endl << htridxc[itd].GetTitle() << endl;
	cout << "  Inside  " << ne << " (" << ne/nb << " per bin)" << endl;
	cout << "  Maximum " << nm << " (factor " << nm/ne*nb << " above mean)" << endl;
	cout << "  at " << xm << endl;
	cout << "  fwhm " << fwhmx << endl;

	TF1 * fgp0x = new TF1( "fgp0x", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", xm-fwhmx, xm+fwhmx );
	fgp0x->SetParameter( 0, nm ); // amplitude
	fgp0x->SetParameter( 1, xm );
	fgp0x->SetParameter( 2, 0.5*fwhmx ); // sigma
	fgp0x->SetParameter( 3, htridxc[itd].GetBinContent(1) ); // BG
	htridxc[itd].Fit( "fgp0x", "qr" );
	cout << "Fit Gauss + BG:"
	     << endl << "  A " << fgp0x->GetParameter(0)
	     << endl << "mid " << fgp0x->GetParameter(1) << " um"
	     << endl << "sig " << fgp0x->GetParameter(2) << " um"
	     << endl << " BG " << fgp0x->GetParameter(3)
	     << endl;

	alignx[ib] -= 0.5*fgp0x->GetParameter(1)*1E-3; // um -> mm
	alignx[ie] -= 0.5*fgp0x->GetParameter(1)*1E-3;

	// dy:

	nb = htridyc[itd].GetNbinsX();
	ne = htridyc[itd].GetSumOfWeights();
	nm = htridyc[itd].GetMaximum();
	double fwhmy = GetFWHM( htridyc[itd] );
	double ym = htridyc[itd].GetBinCenter( htridyc[itd].GetMaximumBin() );
	cout << endl << htridyc[itd].GetTitle() << endl;
	cout << "  Inside  " << ne << " (" << ne/nb << " per bin)" << endl;
	cout << "  Maximum " << nm << " (factor " << nm/ne*nb << " above mean)" << endl;
	cout << "  at " << ym << endl;
	cout << "  fwhm " << fwhmy << endl;

	TF1 * fgp0y = new TF1( "fgp0y", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", ym-fwhmy, ym+fwhmy );
	fgp0y->SetParameter( 0, nm ); // amplitude
	fgp0y->SetParameter( 1, ym );
	fgp0y->SetParameter( 2, 0.5*fwhmy ); // sigma
	fgp0y->SetParameter( 3, htridyc[itd].GetBinContent(1) ); // BG
	htridyc[itd].Fit( "fgp0y", "qr" );
	cout << "Fit Gauss + BG:"
	     << endl << "  A " << fgp0y->GetParameter(0)
	     << endl << "mid " << fgp0y->GetParameter(1) << " um"
	     << endl << "sig " << fgp0y->GetParameter(2) << " um"
	     << endl << " BG " << fgp0y->GetParameter(3)
	     << endl;

	// update:

	aligny[ib] -= 0.5*fgp0y->GetParameter(1)*1E-3; // um -> mm
	aligny[ie] -= 0.5*fgp0y->GetParameter(1)*1E-3;

	// x-y rotation:

	tridxvsy[itd].Fit( "pol1", "q", "", -midy[im], midy[im] );
	TF1 * fdxvsy = tridxvsy[itd].GetFunction( "pol1" );
	cout << endl << tridxvsy[itd].GetTitle()
	     << " slope " << fdxvsy->GetParameter(1)
	     << " = mrad rot"
	     << endl;
	rotx[ib] -= 0.5*fdxvsy->GetParameter(1)*1E-3;
	rotx[ie] -= 0.5*fdxvsy->GetParameter(1)*1E-3;

	tridyvsx[itd].Fit( "pol1", "q", "", -midx[im], midx[im] );
	TF1 * fdyvsx = tridyvsx[itd].GetFunction( "pol1" );
	cout << endl << tridyvsx[itd].GetTitle()
	     << " slope " << fdyvsx->GetParameter(1)
	     << " = mrad rot"
	     << endl;
	roty[ib] += 0.5*fdyvsx->GetParameter(1)*1E-3;
	roty[ie] += 0.5*fdyvsx->GetParameter(1)*1E-3;

      } // itd

    } // aligniteration

    // z-shift:

    if( aligniteration >= 4 ) {

      for( int itd = 0; itd < 2; ++itd ) { // upstream and downstream

	cout << endl;

	int ipl = 3+3*itd; // 3 or 6
	/*
	  if( tridxCvsax[itd].GetEntries() > 999 ) {

	  tridxCvsax[itd].Fit( "pol1", "q", "", -1, 1 ); // um vs mrad
	  TF1 * f1 = tridxCvsax[itd].GetFunction( "pol1" );
	  alignz[ipl-2] += 0.5*f1->GetParameter(1);
	  alignz[ipl]   += 0.5*f1->GetParameter(1);
	  cout << tridxCvsax[itd].GetTitle()
	  << " dz " << f1->GetParameter(1)
	  << ", plane " << ipl-2
	  << " new zpos " << zz[ipl-2] + alignz[ipl-2]
	  << ", plane " << ipl
	  << " new zpos " << zz[ipl] + alignz[ipl]
	  << endl;

	  }
	*/
	if( tridxvstx[itd].GetEntries() > 999 ) {

	  tridxvstx[itd].Fit( "pol1", "q", "", -1, 1 ); // um vs mrad
	  TF1 * f1 = tridxvstx[itd].GetFunction( "pol1" );
	  alignz[ipl-2] -= 0.5*f1->GetParameter(1);
	  alignz[ipl]   -= 0.5*f1->GetParameter(1);
	  cout << tridxvstx[itd].GetTitle()
	       << " dz " << f1->GetParameter(1)
	       << ", plane " << ipl-2
	       << " new zpos " << zz[ipl-2] + alignz[ipl-2]
	       << ", plane " << ipl
	       << " new zpos " << zz[ipl] + alignz[ipl]
	       << endl;

	}

      } // itd

    } // aligniteration

    // driplet vs triplet:

    if( aligniteration >= 5 && hsixdx.GetMaximum() > 999 ) {

      double nb = hsixdx.GetNbinsX();
      double ne = hsixdx.GetSumOfWeights();
      double nm = hsixdx.GetMaximum();
      double fwhmx = GetFWHM( hsixdx );
      double xm = hsixdx.GetBinCenter( hsixdx.GetMaximumBin() );

      cout << endl << hsixdx.GetTitle() << endl;
      cout << "  Inside  " << ne << " (" << ne/nb << " per bin)" << endl;
      cout << "  Maximum " << nm << " (factor " << nm/ne*nb << " above mean)" << endl;
      cout << "  at " << xm << endl;
      cout << "  fwhm " << fwhmx << endl;

      TF1 * fgp0x = new TF1( "fgp0x", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", xm-fwhmx, xm+fwhmx );
      fgp0x->SetParameter( 0, nm ); // amplitude
      fgp0x->SetParameter( 1, xm );
      fgp0x->SetParameter( 2, 0.5*fwhmx ); // sigma
      fgp0x->SetParameter( 3, hsixdx.GetBinContent(1) ); // BG
      hsixdx.Fit( "fgp0x", "qr" );
      cout << "Fit Gauss + BG:"
	   << endl << "  A " << fgp0x->GetParameter(0)
	   << endl << "mid " << fgp0x->GetParameter(1)*1E3 << " um"
	   << endl << "sig " << fgp0x->GetParameter(2)*1E3 << " um"
	   << endl << " BG " << fgp0x->GetParameter(3)
	   << endl;

      // dy:

      nb = hsixdy.GetNbinsX();
      ne = hsixdy.GetSumOfWeights();
      nm = hsixdy.GetMaximum();
      double fwhmy = GetFWHM( hsixdy );
      double ym = hsixdy.GetBinCenter( hsixdy.GetMaximumBin() );

      cout << endl << hsixdy.GetTitle() << endl;
      cout << "  Inside  " << ne << " (" << ne/nb << " per bin)" << endl;
      cout << "  Maximum " << nm << " (factor " << nm/ne*nb << " above mean)" << endl;
      cout << "  at " << ym << endl;
      cout << "  fwhm " << fwhmy << endl;

      TF1 * fgp0y = new TF1( "fgp0y", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", ym-fwhmy, ym+fwhmy );
      fgp0y->SetParameter( 0, nm ); // amplitude
      fgp0y->SetParameter( 1, ym );
      fgp0y->SetParameter( 2, 0.5*fwhmy ); // sigma
      fgp0y->SetParameter( 3, hsixdy.GetBinContent(1) ); // BG
      hsixdy.Fit( "fgp0y", "qr" );
      cout << "Fit Gauss + BG:"
	   << endl << "  A " << fgp0y->GetParameter(0)
	   << endl << "mid " << fgp0y->GetParameter(1)*1E3 << " um"
	   << endl << "sig " << fgp0y->GetParameter(2)*1E3 << " um"
	   << endl << " BG " << fgp0y->GetParameter(3)
	   << endl;

      // update driplet planes:

      for( int ipl = 4; ipl <= 6; ++ipl ) {
	alignx[ipl] += fgp0x->GetParameter(1); // [mm]
	aligny[ipl] += fgp0y->GetParameter(1);
      }

      delete fgp0x;
      delete fgp0y;

    } // aligniteration

    if( aligniteration >= 6 && hsixdx.GetMaximum() > 999 ) {

      // x-y rotation from profiles:

      if( sixdxvsy.GetEntries() > 999 ) {
	sixdxvsy.Fit( "pol1", "q", "", -midy[2], midy[2] );
	TF1 * fdxvsy = sixdxvsy.GetFunction( "pol1" );
	cout << endl << sixdxvsy.GetTitle()
	     << " slope " << fdxvsy->GetParameter(1)
	     << " mu/mm rotx" << endl;
	for( int ipl = 4; ipl <= 6; ++ipl )
	  rotx[ipl] += fdxvsy->GetParameter(1)*1E-3; // um -> mm
      }

      if( sixdyvsx.GetEntries() > 999 ) {
	sixdyvsx.Fit( "pol1", "q", "", -midx[2], midx[2] );
	TF1 * fdyvsx = sixdyvsx.GetFunction( "pol1" );
	cout << endl << sixdyvsx.GetTitle()
	     << " slope " << fdyvsx->GetParameter(1)
	     << " mu/mm roty" << endl;
	for( int ipl = 4; ipl <= 6; ++ipl )
	  roty[ipl] -= fdyvsx->GetParameter(1)*1E-3; // sign
      }

    } // aligniteration

    // dz from dx vs tx:

    if( aligniteration >= 7 && sixdxvstx.GetEntries() > 999 ) {

      sixdxvstx.Fit( "pol1", "q", "", -1, 1 ); // [mu/mrad]
      TF1 * f1 = sixdxvstx.GetFunction( "pol1" );
      cout << endl << sixdxvstx.GetTitle()
	   << "  dz " << f1->GetParameter(1)
	   << endl;

      for( int ipl = 4; ipl <= 6; ++ipl )
	alignz[ipl] += f1->GetParameter(1);

    } // aligniteration

  } // aligniterations

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // write alignment to file:

  if( pbeam > 2.7 ) {

    ofstream alignFile( alignFileName.str() );

    alignFile << "# telescope alignment for run " << run << endl;

    alignFile << "iteration " << aligniteration << endl;

    for( int ipl = 1; ipl <= 6; ++ipl ) {
      alignFile << endl;
      alignFile << "plane " << ipl << endl;
      alignFile << "shiftx " << alignx[ipl] << endl;
      alignFile << "shifty " << aligny[ipl] << endl;
      alignFile << "shiftz " << alignz[ipl] << endl;
      alignFile << "rotxvsy " << rotx[ipl] << endl;
      alignFile << "rotyvsx " << roty[ipl] << endl;
    } // ipl

    alignFile.close();

    cout << endl
	 << "wrote telescope alignment iteration " << aligniteration
	 << " to " << alignFileName.str()
	 << endl;
    if( aligniteration <= 8 )
      cout << endl << "need more align iterations: please run again!" << endl;

  } // high p

  cout << endl << histoFile->GetName() << endl;

  //delete histoFile; // does not help against final crash

  cout << endl;

  return 0;
}
