
// Daniel Pitzl, DESY, Jun 2018, Apr 2019
// telescope analysis with RD53A: edge-on
// module in front

// make edg53
// edg53 36094
// needs runs.dat
// needs geo_201x_yy.dat
// needs align_36094.dat from tele
// uses hot_36094.dat from tele
// uses alignDUT_36094.dat
// uses alignMOD_36094.dat

#include "eudaq/FileReader.hh"
#include "eudaq/PluginManager.hh"
#include "../main/lib/plugins/BDAQ53ConverterPlugin.cc"

#include <TFile.h>
#include <TH1I.h> // counting
#include <TH1D.h> // weighted counts
#include <TH2I.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TF1.h>

#include <sstream> // stringstream
#include <fstream> // filestream
#include <set>
#include <cmath>

using namespace std;
using namespace eudaq;

struct pixel {
  int col;
  int row;
  int tot;
  int frm;
};

struct cluster {
  vector <pixel> vpix; // Armin Burgmeier: list
  int size;
  int ncol, nrow, nfrm;
  double col, row;
  int signal;
  double mindxy;
};

struct triplet {
  double xm;
  double ym;
  double zm;
  double sx;
  double sy;
  bool lk;
  double ttdmin;
  unsigned iA;
  unsigned iB;
  unsigned iC;
};

//------------------------------------------------------------------------------
vector <cluster> getClusn( vector <pixel> pb, int fCluCut = 1 ) // 1 = no gap
{
  // returns clusters with pixel coordinates
  // next-neighbour topological clustering (allows fCluCut-1 empty pixels)

  vector <cluster> vc;
  if( pb.size() == 0 ) return vc;

  int* gone = new int[pb.size()];
  for( unsigned i = 0; i < pb.size(); ++i )
    gone[i] = 0;

  unsigned seed = 0;

  while( seed < pb.size() ) {

    // start a new cluster

    cluster c;
    c.vpix.push_back( pb[seed] );
    gone[seed] = 1;

    // let it grow as much as possible:

    int growing;
    do {
      growing = 0;
      for( unsigned i = 0; i < pb.size(); ++i ) {
        if( !gone[i] ){ // unused pixel
          for( unsigned int p = 0; p < c.vpix.size(); ++p ) { // vpix in cluster so far
            int dr = c.vpix.at(p).row - pb[i].row;
            int dc = c.vpix.at(p).col - pb[i].col;
            if( (   dr>=-fCluCut) && (dr<=fCluCut) 
		&& (dc>=-fCluCut) && (dc<=fCluCut) ) {
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
	  ++p->tot;
	  ++q->tot;
	}
    }

    // added all I could. determine position and append it to the list of clusters:

    c.size = c.vpix.size();
    c.col = 0;
    c.row = 0;
    double sumnn = 0;
    int minx = 999;
    int maxx = 0;
    int miny = 999;
    int maxy = 0;
    int minf = 999;
    int maxf = 0;

    for( vector<pixel>::iterator p = c.vpix.begin(); p != c.vpix.end(); ++p ) {

      int nn = max( 1, p->tot ); // neighbours
      sumnn += nn;
      c.col += p->col * nn;
      c.row += p->row * nn;

      if( p->col > maxx ) maxx = p->col;
      if( p->col < minx ) minx = p->col;
      if( p->row > maxy ) maxy = p->row;
      if( p->row < miny ) miny = p->row;
      if( p->frm > maxf ) maxf = p->frm;
      if( p->frm < minf ) minf = p->frm;

    }

    //cout << "(cluster with " << c.vpix.size() << " pixels)" << endl;

    c.col /= sumnn;
    c.row /= sumnn;
    c.signal = sumnn;
    c.ncol = maxx-minx+1;
    c.nrow = maxy-miny+1;
    c.nfrm = maxf-minf+1;
    c.mindxy = 999;

    vc.push_back(c); // add cluster to vector

    // look for a new seed = used pixel:

    while( ( ++seed < pb.size() ) && gone[seed] );

  } // while over seeds

  delete[] gone;

  return vc; // vector of clusters

} // getclusn

//------------------------------------------------------------------------------
vector <cluster> getClusq( vector <pixel> pb, int fCluCut = 1 ) // 1 = no gap
{
  // returns clusters with pixel coordinates
  // next-neighbour topological clustering (allows fCluCut-1 empty pixels)

  vector <cluster> vc;
  if( pb.size() == 0 ) return vc;

  int* gone = new int[pb.size()];
  for( unsigned i = 0; i < pb.size(); ++i )
    gone[i] = 0;

  unsigned seed = 0;

  while( seed < pb.size() ) {

    // start a new cluster

    cluster c;
    c.vpix.push_back( pb[seed] );
    gone[seed] = 1;

    // let it grow as much as possible:

    int growing;
    do {
      growing = 0;
      for( unsigned i = 0; i < pb.size(); ++i ) {
        if( !gone[i] ){ // unused pixel
          for( unsigned int p = 0; p < c.vpix.size(); ++p ) { // vpix in cluster so far
            int dr = c.vpix.at(p).row - pb[i].row;
            int dc = c.vpix.at(p).col - pb[i].col;
            if( (   dr>=-fCluCut) && (dr<=fCluCut) 
		&& (dc>=-fCluCut) && (dc<=fCluCut) ) {
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

    // added all I could. determine position and append it to the list of clusters:

    c.size = c.vpix.size();
    c.col = 0;
    c.row = 0;
    double sumQ = 0;
    int minx = 999;
    int maxx = 0;
    int miny = 999;
    int maxy = 0;
    int minf = 999;
    int maxf = 0;

    for( vector<pixel>::iterator p = c.vpix.begin(); p != c.vpix.end(); ++p ) {

      double Qpix = p->tot;

      sumQ += Qpix;

      c.col += p->col*Qpix;
      c.row += p->row*Qpix;

      if( p->col > maxx ) maxx = p->col;
      if( p->col < minx ) minx = p->col;
      if( p->row > maxy ) maxy = p->row;
      if( p->row < miny ) miny = p->row;
      if( p->frm > maxf ) maxf = p->frm;
      if( p->frm < minf ) minf = p->frm;

    }

    //cout << "(cluster with " << c.vpix.size() << " pixels)" << endl;

    if( sumQ > 0 ) {
      c.col /= sumQ;
      c.row /= sumQ;
    }
    else {
      c.col = (*c.vpix.begin()).col;
      c.row = (*c.vpix.begin()).row;
      cout << "GetClus: cluster with signal" << sumQ << endl;
    }

    c.signal = sumQ;
    c.ncol = maxx-minx+1;
    c.nrow = maxy-miny+1;
    c.nfrm = maxf-minf+1;
    c.mindxy = 999;

    vc.push_back(c); // add cluster to vector

    // look for a new seed = used pixel:

    while( ( ++seed < pb.size() ) && gone[seed] );

  } // while over seeds

  delete gone;

  return vc; // vector of clusters

} // getclusq

//------------------------------------------------------------------------------
int main( int argc, char* argv[] )
{
  cout << "main " << argv[0] << " called with " << argc << " arguments" << endl;

  if( argc == 1 ) {
    cout << "give run number" << endl;
    return 1;
  }

  // run number = last arg

  string runnum( argv[argc-1] );
  int run = atoi( argv[argc-1] );

  cout << "run " << run << endl;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // further arguments:

  int fev = 0; // 1st event
  int lev = 999222111; // last event

  for( int i = 1; i < argc; ++i ) {

    if( !strcmp( argv[i], "-f" ) )
      fev = atoi( argv[++i] ); // 1st event

    if( !strcmp( argv[i], "-l" ) )
      lev = atoi( argv[++i] ); // last event

  } // argc

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // runs.dat:

  cout << endl;

  string geoFileName( "geo.dat" );
  double pbeam = 5.2;
  int chip0 = 501;
  int fifty = 0; // default is 100x25
  int modrun = 0;

  ifstream runsFile( "runs.dat" );

  if( runsFile.bad() || ! runsFile.is_open() ) {
    cout << "Error opening runs.dat" << endl;
    return 1;
  }
  // can there be instructions between if and else ? no

  else {

    cout << "read runs.dat:" << endl;

    string hash( "#" );
    string RUN( "run" );
    string MODRUN( "modrun" );
    string GEO( "geo" );
    string GeV( "GeV" );
    string CHIP( "chip" );
    string TURN( "turn" );
    string TILT( "tilt" );
    string FIFTY( "fifty" );
    bool found = 0;

    while( ! runsFile.eof() ) {

      string line;
      getline( runsFile, line );

      if( line.empty() ) continue;

      istringstream tokenizer( line );
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

      if( tag == MODRUN ) {
	tokenizer >> modrun;
	continue;
      }

      if( tag == GEO ) {
	tokenizer >> geoFileName;
	continue;
      }

      if( tag == GeV ) {
	tokenizer >> pbeam;
	continue;
      }

      if( tag == CHIP ) {
	tokenizer >> chip0;
	continue;
      }

      if( tag == FIFTY ) {
	tokenizer >> fifty;
	continue;
      }

      // anything else on the line and in the file gets ignored

    } // while getline

    if( found )
      cout 
	<< "  beam " << pbeam << " GeV" << endl
	<< "  geo file " << geoFileName << endl
	<< "  DUT chip " << chip0 << endl
	<< "  fifty " << fifty << endl
	<< "  modrun " << modrun << endl
	;
    else {
      cout << "run " << run << " not found in runs.dat" << endl;
      return 1;
    }

  } // runsFile

  runsFile.close();

  const double fTLU = 384E6; // 384 MHz TLU clock

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

  cout << endl;

  if( geoFile.bad() || ! geoFile.is_open() ) {
    cout << "Error opening " << geoFileName << endl;
    return 1;
  }

  cout << "read geometry from " << geoFileName << endl;

  { // open local scope

    string hash( "#" );
    string plane( "plane" );
    string type( "type" );
    string sizexs( "sizex" );
    string sizeys( "sizey" );
    string npixelx( "npixelx" );
    string npixely( "npixely" );
    string zpos( "zpos" );

    int ipl = 0;
    string chiptype;

    while( ! geoFile.eof() ) {

      string line;
      getline( geoFile, line );
      cout << line << endl;

      if( line.empty() ) continue;

      istringstream tokenizer( line );
      string tag;
      tokenizer >> tag; // leading white space is suppressed
      if( tag.substr(0,1) == hash ) // comments start with #
	continue;

      if( tag == plane ) {
	tokenizer >> ipl;
	continue;
      }

      if( ipl < 0 || ipl > 8 ) {
	cout << "geo wrong plane number " << ipl << endl;
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
      ptchx[ipl] = sizex[ipl] / nx[ipl]; // pixel size
      ptchy[ipl] = sizey[ipl] / ny[ipl];
      midx[ipl] = 0.5 * sizex[ipl]; // mid of plane
      midy[ipl] = 0.5 * sizey[ipl]; // mid of plane
    }

  } // geo scope

  cout << endl;
  for( int ipl = 1; ipl <= 6; ++ipl )
    cout << ipl << " zz " << zz[ipl] << endl;

  geoFile.close();

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // read Mimosa telescope alignment:

  int aligniteration = 0;
  double alignx[9];
  double aligny[9];
  double alignz[9];
  double rotx[9];
  double roty[9];

  ostringstream alignFileName; // output string stream

  alignFileName << "align_" << run << ".dat";

  ifstream ialignFile( alignFileName.str() );

  cout << endl;

  if( ialignFile.bad() || ! ialignFile.is_open() ) {
    cout << "Error opening " << alignFileName.str() << endl
	 << "  please do: tele -g " << geoFileName << " " << run << endl
	 << endl;
    return 1;
  }
  else {

    cout << "read alignment from " << alignFileName.str() << endl;

    string hash( "#" );
    string iteration( "iteration" );
    string plane( "plane" );
    string shiftx( "shiftx" );
    string shifty( "shifty" );
    string shiftz( "shiftz" );
    string rotxvsy( "rotxvsy" );
    string rotyvsx( "rotyvsx" );

    int ipl = 1;

    while( ! ialignFile.eof() ) {

      string line;
      getline( ialignFile, line );
      cout << line << endl;

      if( line.empty() ) continue;

      istringstream tokenizer( line );
      string tag;
      tokenizer >> tag; // leading white space is suppressed
      if( tag.substr(0,1) == hash ) // comments start with #
	continue;

      if( tag == iteration ) 
	tokenizer >> aligniteration;

      if( tag == plane )
	tokenizer >> ipl;

      if( ipl < 1 || ipl > 6 ) { // Mimosa
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

  cout << endl;
  for( int ipl = 1; ipl <= 6; ++ipl )
    cout << ipl << " alignz " << alignz[ipl] << endl;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // (re-)create root file:

  ostringstream rootFileName; // output string stream

  rootFileName << "edg53_" << run << ".root";

  TFile* histoFile = new TFile( rootFileName.str(  ).c_str(  ), "RECREATE" );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // telescope hot pixels:

  ostringstream hotFileName; // output string stream

  hotFileName << "hot_" << run << ".dat";

  ifstream ihotFile( hotFileName.str() );

  set <int> hotset[9];

  cout << endl;

  if( ihotFile.bad() || ! ihotFile.is_open() ) {
    cout << "no " << hotFileName.str() << " (created by tele)" << endl;
  }
  else {

    cout << "read hot pixel list from " << hotFileName.str() << endl;

    string hash( "#" );
    string plane( "plane" );
    string pix( "pix" );

    int ipl = 0;

    while( ! ihotFile.eof() ) {

      string line;
      getline( ihotFile, line );
      //cout << line << endl;

      if( line.empty() ) continue;

      istringstream tokenizer( line );
      string tag;
      tokenizer >> tag; // leading white space is suppressed
      if( tag.substr(0,1) == hash ) // comments start with #
	continue;

      if( tag == plane )
	tokenizer >> ipl;

      if( ipl < 1 || ipl > 6 ) { // Mimosa
	cout << "hot wrong plane number " << ipl << endl;
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

  for( int ipl = 0; ipl <= 6; ++ipl )
    cout << "  plane " << ipl << ": hot " << hotset[ipl].size() << endl;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // DUT:

  const double wt = atan(1.0) / 45.0; // pi/180 deg

  int iDUT = 0; // eudaq

  int DUTaligniteration = 0;
  double DUTalignx = 0.0;
  double DUTaligny = 0.0;
  double DUTrot = 0; // [rad]
  double DUTtilt = 0; // [rad]
  double DUTturn = 0; // [rad]
  double DUTz = 0.5 * ( zz[3] + zz[4] );

  ostringstream DUTalignFileName; // output string stream

  DUTalignFileName << "alignDUT_" << run << ".dat";

  ifstream iDUTalignFile( DUTalignFileName.str() );

  cout << endl;

  if( iDUTalignFile.bad() || ! iDUTalignFile.is_open() ) {
    cout << "no " << DUTalignFileName.str() << ", will bootstrap" << endl;
  }
  else {

    cout << "read DUTalignment from " << DUTalignFileName.str() << endl;

    string hash( "#" );
    string iteration( "iteration" );
    string alignx( "alignx" );
    string aligny( "aligny" );
    string rot( "rot" );
    string tilt( "tilt" );
    string turn( "turn" );
    string dz( "dz" );

    while( ! iDUTalignFile.eof() ) {

      string line;
      getline( iDUTalignFile, line );
      cout << line << endl;

      if( line.empty() ) continue;

      istringstream tokenizer( line );
      string tag;
      tokenizer >> tag; // leading white space is suppressed
      if( tag.substr(0,1) == hash ) // comments start with #
	continue;

      if( tag == iteration ) 
	tokenizer >> DUTaligniteration;

      double val;
      tokenizer >> val;
      if(      tag == alignx )
	DUTalignx = val;
      else if( tag == aligny )
	DUTaligny = val;
      else if( tag == rot )
	DUTrot = val; // [rad]
      else if( tag == tilt )
	DUTtilt = val; // [rad]
      else if( tag == turn )
	DUTturn = val; // [rad]
      else if( tag == dz )
	DUTz = val + zz[3];

      // anything else on the line and in the file gets ignored

    } // while getline

  } // alignFile

  iDUTalignFile.close();

  double cf = cos( DUTrot );
  double sf = sin( DUTrot );
  double ca = cos( DUTtilt );
  double sa = sin( DUTtilt );
  double co = cos( DUTturn );
  double so = sin( DUTturn );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // DUT hot pixels:

  cout << endl;

  ostringstream DUThotFileName; // output string stream

  DUThotFileName << "hotDUT_" << run << ".dat";

  ifstream iDUThotFile( DUThotFileName.str() );

  if( iDUThotFile.bad() || ! iDUThotFile.is_open() ) {
    cout << "no " << DUThotFileName.str() << endl;
  }
  else {

    cout << "read DUT hot pixel list from " << DUThotFileName.str() << endl;

    string hash( "#" );
    string pix( "pix" );

    while( ! iDUThotFile.eof() ) {

      string line;
      getline( iDUThotFile, line );
      //cout << line << endl;

      if( line.empty() ) continue;

      stringstream tokenizer( line );
      string tag;
      tokenizer >> tag; // leading white space is suppressed
      if( tag.substr(0,1) == hash ) // comments start with #
	continue;

      if( tag == pix ) {

	int ix, iy;
	tokenizer >> ix; // ROC col
	tokenizer >> iy; // ROC row
	int ipx = ix * ny[iDUT] + iy;
	hotset[iDUT].insert(ipx);

      }

    } // while getline

  } // hotFile

  iDUThotFile.close();

  cout << "DUT hot " << hotset[iDUT].size() << endl;

  cout << endl;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // MOD:

  string modfileA = "mod/run" + to_string(modrun) + "A.out"; // C++11
  cout << "try to open  " << modfileA;
  ifstream Astream( modfileA.c_str() );
  if( !Astream ) {
    cout << " : failed " << endl;
    modrun = 0; // flag
  }
  else
    cout << " : succeed " << endl;

  string modfileB = "mod/run" + to_string(modrun) + "B.out";
  cout << "try to open  " << modfileB;
  ifstream Bstream( modfileB.c_str() );
  if( !Bstream ) {
    cout << " : failed " << endl;
    modrun = 0;
  }
  else
    cout << " : succeed " << endl;

  string sl;
  int mev = 0;
  if( fev ) cout << "MOD skip " << fev << endl;
  while( mev < fev ) {
    getline( Astream, sl ); // fast forward
    getline( Bstream, sl );
    ++mev;
  }

  int iMOD = 7;

  int MODaligniteration = 0;
  double MODalignx = 0;
  double MODaligny = 0;
  double MODrot = 0;
  double MODtilt =  19.2; // [deg]
  double MODturn = -27.0; // [deg]
  double MODz = -45 + zz[1];

  ostringstream MODalignFileName; // output string stream

  MODalignFileName << "alignMOD_" << run << ".dat";

  ifstream iMODalignFile( MODalignFileName.str() );

  cout << endl;

  if( iMODalignFile.bad() || ! iMODalignFile.is_open() ) {
    cout << "no " << MODalignFileName.str() << ", will bootstrap" << endl;
  }
  else {

    cout << "read MODalignment from " << MODalignFileName.str() << endl;

    string hash( "#" );
    string iteration( "iteration" );
    string alignx( "alignx" );
    string aligny( "aligny" );
    string rot( "rot" );
    string tilt( "tilt" );
    string turn( "turn" );
    string dz( "dz" );

    while( ! iMODalignFile.eof() ) {

      string line;
      getline( iMODalignFile, line );
      cout << line << endl;

      if( line.empty() ) continue;

      istringstream tokenizer( line );
      string tag;
      tokenizer >> tag; // leading white space is suppressed
      if( tag.substr(0,1) == hash ) // comments start with #
	continue;

      if( tag == iteration ) 
	tokenizer >> MODaligniteration;

      double val;
      tokenizer >> val;
      if(      tag == alignx )
	MODalignx = val;
      else if( tag == aligny )
	MODaligny = val;
      else if( tag == rot )
	MODrot = val;
      else if( tag == tilt )
	MODtilt = val;
      else if( tag == turn )
	MODturn = val;
      else if( tag == dz )
	MODz = val + zz[1];

      // anything else on the line and in the file gets ignored

    } // while getline

  } // alignFile

  iMODalignFile.close();

  // normal vector on MOD surface:
  // N = ( 0, 0, -1 ) on MOD, towards -z
  // transform into tele system:
  // tilt alpha around x
  // turn omega around y

  const double com = cos( MODturn*wt );
  const double som = sin( MODturn*wt );
  const double cam = cos( MODtilt*wt );
  const double sam = sin( MODtilt*wt );
  const double cfm = cos( MODrot );
  const double sfm = sin( MODrot );

  const double Nxm =-cam*som;
  const double Nym = sam;
  const double Nzm =-cam*com;

  const double normm = cos( MODturn*wt ) * cos( MODtilt*wt ); // length of Nz

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // book histos:

  double f = 5.2 / pbeam;

  TH1I t1Histo( "t1", "event time;event time [s];events / 10 ms", 100, 0, 1 );
  TH1I t2Histo( "t2", "event time;event time [s];events / s", 300, 0, 300 );
  TH1I t3Histo( "t3", "event time;event time [s];events / 10 s", 150, 0, 1500 );
  TH1I t4Histo( "t4", "event time;event time [s];events /10 s", 600, 0, 6000 );
  TH1I t5Histo( "t5", "event time;event time [s];events / 100 s", 600, 0, 60000 );
  TH1I t6Histo( "t6", "event time;event time [h];events / 3 min", 1000, 0, 50 );

  TH1I dtusHisto( "dtus", "time between events;time between events [us];events", 100, 0, 1000 );
  TH1I dtmsHisto( "dtms", "time between events;time between events [ms];events", 100, 0, 1000 );
  TH1I dt373Histo( "dt373", "time between events;time between events mod 373;events", 400, -0.5, 399.5 );
  TH1I dt374Histo( "dt374", "time between events;time between events mod 374;events", 400, -0.5, 399.5 );
  TH1I dt375Histo( "dt375", "time between events;time between events mod 375;events", 400, -0.5, 399.5 );
  TH1I dt376Histo( "dt376", "time between events;time between events mod 376;events", 400, -0.5, 399.5 );
  TH1I dt377Histo( "dt377", "time between events;time between events mod 377;events", 400, -0.5, 399.5 );
  TProfile dt375vsdt( "dt375vsdt", "dt vs dt;time between events [us];<time between events mod 375>",
		      200, 0, 2000, 0, 400 );

  TH1I hpivot[9];
  TH1I hnpx[9];
  TH1I hnframes[9];
  TH1I hcol[9];
  TH1I hrow[9];
  TH2I * hmap[9];

  TH1I hncl[9];
  TH1I hsiz[9];
  TH1I hncol[9];
  TH1I hnrow[9];
  TH1I hdxy[9];

  for( int ipl = 0; ipl <= 7; ++ipl ) {

    int nbx = 400; // RD53
    int nby = 192;
    int mx = nbx;
    int my = nby;

    if( ipl >= 1 ) {
      nbx = nx[1]/2; // Mimosa
      nby = ny[1]/2; // for next round
      mx = nx[1];
      my = ny[1];
    }
    if( ipl == 7 ) { // MOD
      nbx = 8*(52+2);
      nby = 162; // with big pix
      mx = nbx;
      my = nby;
    }

    hpivot[ipl] = TH1I( Form( "pivot%i", ipl ),
			Form( "%i pivot;pivot;plane %i events", ipl, ipl ),
			1000, 0, 10*1000 );

    hnpx[ipl] = TH1I( Form( "npx%i", ipl ),
		      Form( "%i pixel per event;pixels;plane %i events", ipl, ipl ),
		      201, -0.5, 200.5 );

    hnframes[ipl] = TH1I( Form( "nframes%i", ipl ),
			  Form( "%i frames per event;frames;plane %i events", ipl, ipl ),
			  51, -0.5, 50.5 );

    hcol[ipl] = TH1I( Form( "col%i", ipl ),
		      Form( "%i col;col;plane %i pixels", ipl, ipl ), 
		      nbx, 0, mx );
    hrow[ipl] = TH1I( Form( "row%i", ipl ),
		      Form( "%i row;row;plane %i pixels", ipl, ipl ),
		      nby, 0, my );
    hmap[ipl] = new TH2I( Form( "map%i", ipl ),
			  Form( "%i map;col;row;plane %i pixels", ipl, ipl ),
			  nbx, 0, mx, nby, 0, my );

    hncl[ipl] = TH1I( Form( "ncl%i", ipl ),
		      Form( "plane %i cluster per event;cluster;plane %i events", ipl, ipl ),
		      51, -0.5, 50.5 );
    hsiz[ipl] = TH1I( Form( "clsz%i", ipl ),
		      Form( "%i cluster size;pixels/cluster;plane %i clusters", ipl, ipl ),
		      51, -0.5, 50.5 );
    hncol[ipl] = TH1I( Form( "ncol%i", ipl ), 
		       Form( "%i cluster size x;columns/cluster;plane %i clusters", ipl, ipl ),
		       21, -0.5, 20.5 );
    hnrow[ipl] = TH1I( Form( "nrow%i", ipl ),
		       Form( "%i cluster size y;rows/cluster;plane %i clusters", ipl, ipl ),
		       21, -0.5, 20.5 );
    hdxy[ipl] = TH1I( Form( "dxy%i", ipl ),
		      Form( "%i cluster distance;cluster distance [pixels];plane %i clusters", ipl, ipl ),
		      100, 0, 10 );

  } // planes

  TProfile dutnpxvsev( "dutnpxvsev",
		       "DUT pixels vs events;events;DUT pixels / 1000",
		       500, 0, 500E3 );

  TH1I dutpxbc0Histo( "dutpxbc0",
		    "DUT pixel BC;DUT pixel BC;DUT pixels without masking",
		    32, -0.5, 31.5 );

  hmap[8] = new TH2I( Form( "map%i", 8 ),
		      Form( "%i map;col;row;plane %i pixels", 8, 8 ),
		      400, 0, 400, 192, 0, 192 );

  TH1I dutpxcolHisto( "dutpxcol",
		      "DUT pixel column;DUT pixel column;DUT pixels",
		      nx[iDUT], 0, nx[iDUT] );
  TH1I dutpxrowHisto( "dutpxrow",
		      "DUT pixel row;DUT pixel row;DUT pixels",
		      ny[iDUT], 0, ny[iDUT] );
  TH1I dutpxqHisto( "dutpxq",
		    "DUT pixel signal;DUT pixel signal [ToT];DUT pixels",
		    16, 0, 16 );
  TH1I dutpxbcHisto( "dutpxbc",
		    "DUT pixel BC;DUT pixel BC;DUT pixels",
		    32, 0, 32 );

  TProfile dutpxqvsx( "dutpxqvsx",
		      "DUT pixel signal vs x;column;<pixel signal> [ToT]",
		      400, 0, 400, 0, 16 );
  TProfile2D * dutpxqvsxy = new
    TProfile2D( "dutpxqvsxy",
		"DUT pixel signal map;column;row;<pixel signal> [ToT]",
		400, 0, 400, 192, 0, 192, 0, 16 );

  // triplets:

  TH1I hdx13( "dx13", "1-3 dx;1-3 dx [mm];cluster pairs", 100, -f, f );
  TH1I hdy13( "dy13", "1-3 dy;1-3 dy [mm];cluster pairs", 100, -f, f );

  TH1I htridx( "tridx", "triplet dx;triplet dx [mm];triplets", 100, -0.1, 0.1 );
  TH1I htridy( "tridy", "triplet dy;triplet dy [mm];triplets", 100, -0.1, 0.1 );

  TH1I htridxc( "tridxc", "triplet dx;triplet dx [mm];triplets", 100, -0.05, 0.05 );
  TH1I htridyc( "tridyc", "triplet dy;triplet dy [mm];triplets", 100, -0.05, 0.05 );

  TProfile tridxvsx( "tridxvsx",
		     "triplet dx vs x;triplet x [mm];<triplet #Deltax> [mm]",
		     120, -12, 12, -0.05, 0.05 );
  TProfile tridxvsy( "tridxvsy",
		     "triplet dx vs y;triplet y [mm];<triplet #Deltax> [mm]",
		     110, -5.5, 5.5, -0.05, 0.05 );
  TProfile tridxvstx( "tridxvstx",
		      "triplet dx vs slope x;triplet slope x [rad];<triplet #Deltax> [mm]",
		      60, -0.003, 0.003, -0.05, 0.05 );
  TProfile tridxvst3( "tridxvst3",
		      "triplet dx vs time;time [s];<triplet #Deltax> [mm] / 10s",
		      300, 0, 3000, -0.05, 0.05 );
  TProfile tridxvst5( "tridxvst5",
		      "triplet dx vs time;time [s];<triplet #Deltax> [mm] / min",
		      1100, 0, 66000, -0.05, 0.05 );

  TProfile tridyvsx( "tridyvsx",
		     "triplet dy vs x;triplet x [mm];<triplet #Deltay> [mm]",
		     110, -11, 11, -0.05, 0.05 );
  TProfile tridyvsty( "tridyvsty",
		      "triplet dy vs slope y;triplet slope y [rad];<triplet #Deltay> [mm]",
		      60, -0.003, 0.003, -0.05, 0.05 );
  TProfile tridyvst3( "tridyvst3",
		      "triplet dy vs time;time [s];<triplet #Deltay> [mm] / 10s",
		      300, 0, 3000, -0.05, 0.05 );
  TProfile tridyvst5( "tridyvst5",
		      "triplet dy vs time;time [h];<triplet #Deltay> [mm] / min",
		      1100, 0, 66000, -0.05, 0.05 );

  TH1I trixHisto( "trix", "triplets x;x [mm];triplets",
		  240, -12, 12 );
  TH1I triyHisto( "triy", "triplets y;y [mm];triplets",
		  120, -6, 6 );
  TH2I * trixyHisto = new
    TH2I( "trixy", "triplets x-y;x [mm];y [mm];triplets",
	  240, -12, 12, 120, -6, 6 );
  TH1I tritxHisto( "tritx", "triplet slope x;slope x [rad];triplets",
		   100, -0.005*f, 0.005*f );
  TH1I trityHisto( "trity", "triplet slope y;slope y [rad];triplets",
		   100, -0.005*f, 0.005*f );

  TH1I ntriHisto( "ntri", "triplets;triplets;events", 51, -0.5, 50.5 );

  // MOD vs triplets:

  TH1I ttdminmod1Histo( "ttdminmod1",
		     "triplet isolation at MOD;triplet at MOD min #Delta_{xy} [mm];triplet pairs",
		     100, 0, 1 );
  TH1I ttdminmod2Histo( "ttdminmod2",
		     "triplet isolation at MOD;triplet at MOD min #Delta_{xy} [mm];triplet pairs",
		     150, 0, 15 );

  TH1I modxHisto( "modx",
		  "MOD x;MOD cluster x [mm];MOD clusters",
		  700, -35, 35 );
  TH1I modyHisto( "mody",
		  "MOD y;MOD cluster y [mm];MOD clusters",
		  200, -10, 10 );

  TH1I modsxaHisto( "modsxa",
		    "MOD + triplet x;MOD cluster + triplet #Sigmax [mm];MOD clusters",
		    1280, -32, 32 );
  TH1I moddxaHisto( "moddxa",
		    "MOD - triplet x;MOD cluster - triplet #Deltax [mm];MOD clusters",
		    1280, -32, 32 );

  TH1I modsyaHisto( "modsya",
		    "MOD + triplet y;MOD cluster + triplet #Sigmay [mm];MOD clusters",
		    320, -8, 8 );
  TH1I moddyaHisto( "moddya",
		    "MOD - triplet y;MOD cluster - triplet #Deltay [mm];MOD clusters",
		    320, -8, 8 );

  TH1I moddxHisto( "moddx",
		   "MOD - triplet x;MOD cluster - triplet #Deltax [mm];MOD clusters",
		   500, -2.5, 2.5 );
  TH1I moddxcHisto( "moddxc",
		    "MOD - triplet x;MOD cluster - triplet #Deltax [mm];MOD clusters",
		    200, -0.5, 0.5 );
  TProfile moddxvsx( "moddxvsx",
		     "MOD #Deltax vs x;x track [mm];<cluster - triplet #Deltax> [mm]",
		     216, -32.4, 32.4, -2.5, 2.5 );
  TProfile moddxvsy( "moddxvsy",
		     "MOD #Deltax vs y;y track [mm];<cluster - triplet #Deltax> [mm]",
		     160, -8, 8, -2.5, 2.5 );
  TProfile moddxvstx( "moddxvstx",
		      "MOD #Deltax vs #theta_{x};x track slope [rad];<cluster - triplet #Deltax> [mm]",
		      80, -0.002, 0.002, -2.5, 2.5 );
  TProfile moddxvst5( "moddxvst5",
		      "MOD dx vs time;time [s];<MOD #Deltax> [mm] / min",
		      1100, 0, 66000, -0.05, 0.05 );

  TH1I moddyHisto( "moddy",
		   "MOD - triplet y;MOD cluster - triplet #Deltay [mm];MOD clusters",
		   200, -0.5, 0.5 );
  TH1I moddycHisto( "moddyc",
		    "MOD - triplet y;MOD cluster - triplet #Deltay [mm];MOD clusters",
		    200, -0.5, 0.5 );
  TH1I moddycqHisto( "moddycq",
		     "MOD - triplet y Landau peak;MOD cluster - triplet #Deltay [mm];Landau peak MOD clusters",
		     500, -0.5, 0.5 );
  TProfile moddyvsx( "moddyvsx",
		     "MOD #Deltay vs x;x track [mm];<cluster - triplet #Deltay> [mm]",
		     216, -32.4, 32.4, -0.5, 0.5 );
  TProfile moddyvsy( "moddyvsy",
		     "MOD #Deltay vs y;y track [mm];<cluster - triplet #Deltay> [mm]",
		     160, -8, 8, -0.5, 0.5 );
  TProfile moddyvsty( "moddyvsty",
		      "MOD #Deltay vs #theta_{y};y track slope [rad];<cluster - triplet #Deltay> [mm]",
		      80, -0.002, 0.002, -0.5, 0.5 );
  TProfile moddyvst5( "moddyvst5",
		      "MOD dy vs time;time [h];<MOD #Deltay> [mm] / min",
		      1100, 0, 66000, -0.05, 0.05 );

  TH1I modnpxHisto( "modnpx",
		    "MOD linked clusters;MOD cluster size [pixels];linked MOD cluster",
		    20, 0.5, 20.5 );

  TH1I modqHisto( "modq",
		  "MOD linked clusters;MOD cluster charge [ke];linked MOD cluster",
		  80, 0, 80 );
  TH1I modq0Histo( "modq0",
		   "MOD linked clusters;MOD normal cluster charge [ke];linked MOD cluster",
		   80, 0, 80 );

  TProfile2D * modnpxvsxmym = new
    TProfile2D( "modnpxvsxmym",
		"MOD cluster size vs xmod ymod;x track mod 300 [#mum];y track mod 200 [#mum];MOD <cluster size> [pixels]",
		120, 0, 300, 80, 0, 200, 0, 20 );

  TH1I modlkxBHisto( "modlkxb",
		     "linked triplet at MOD x;triplet x at MOD [mm];linked triplets",
		     216, -32.4, 32.4 );
  TH1I modlkyBHisto( "modlkyb",
		     "linked triplet at MOD y;triplet y at MOD [mm];linked triplets",
		     160, -8, 8 );
  TH1I modlkxHisto( "modlkx",
		    "linked triplet at MOD x;triplet x at MOD [mm];linked triplets",
		    216, -32.4, 32.4 );
  TH1I modlkyHisto( "modlky",
		    "linked triplet at MOD y;triplet y at MOD [mm];linked triplets",
		    160, -8, 8 );

  TH1I modlkcolHisto( "modlkcol",
		      "MOD linked col;MOD linked col;linked MOD cluster",
		      216, 0, 432 );
  TH1I modlkrowHisto( "modlkrow",
		      "MOD linked row;MOD linked row;linked MOD cluster",
		      182, 0, 182 );

  TH1I ttdmin1Histo( "ttdmin1",
		     "telescope triplets isolation;triplet min #Delta_{xy} [mm];triplet pairs",
		     100, 0, 1 );
  TH1I ttdmin2Histo( "ttdmin2",
		     "telescope triplets isolation;triplet min #Delta_{xy} [mm];triplet pairs",
		     150, 0, 15 );

  // DUT vs triplets:

  TH1I trixcHisto( "trixc",
		   "in-time triplet at DUT x;triplet x at DUT [mm];in-time triplets",
		   120, -6, 6 );
  TH1I triycHisto( "triyc",
		   "in-time triplet at DUT y;triplet y at DUT [mm];in-time triplets",
		   120, -6, 6 );

  TH2I * pixzyHisto = new
    TH2I( "pixzy",
	  "DUT hits;z [mm];y [mm];pixels",
	  220, -11, 11, 480, -6, 6 ); // bin = pix

  TH1I pixdxaHisto( "pixdxa",
		    "PIX - Telescope x;pixel - triplet #Deltax [mm];pixels",
		    440, -11, 11 );
  TH1I pixsxaHisto( "pixsxa",
		    "PIX + Telescope x;pixel + triplet #Sigmax [mm];pixels",
		    440, -11, 11 );
  TH1I pixdyaHisto( "pixdya",
		    "PIX - Telescope y;pixel - triplet #Deltay [mm];pixels",
		    240, -6, 6 );
  TH1I pixsyaHisto( "pixsya",
		    "PIX + Telescope y;pixel + triplet #Sigmay [mm];pixels",
		    240, -6, 6 );

  TH1I pixdyHisto( "pixdy",
		   "pixel - Telescope dy;pixel - triplet #Deltaxy [mm];pixels",
		   200, -0.5, 0.5 );

  TH1I pixdxcHisto( "pixdxc",
		    "pixel - Telescope dx;pixel - triplet #Deltaxx [mm];pixels",
		    200, -0.5, 0.5 );
  TProfile pixdxvsz( "pixdxvsz",
		     "#Deltax vs z;z [mm];<pixel - triplet #Deltax> [mm]",
		     200, -10, 10, -0.1, 0.1 );
  TProfile pixdxvsy( "pixdxvsy",
		     "#Deltax vs y;track y [mm];<pixel - triplet #Deltax> [mm]",
		     100, -5, 5, -0.1, 0.1 );
  TProfile pixdxvsev( "pixdxvsev", "PIX - Telescope x vs time;time [events];<#Deltax> [mm]",
		      600, 0, 600*1000, -0.1, 0.1 );

  TH2I * hroadmap = new TH2I( "roadmap",
			      "road map;col;row;active pixels",
			      400, 0, 400, 192, 0, 192 );
  TH1I pixbclkHisto( "pixbclk",
		     "DUT linked pixel BC;DUT pixel BC;DUT pixels on tracks",
		     32, 0, 32 );
  TProfile pixbcvsd( "pixbcvsd",
		     "PIX BC vs depth;depth [mm];<pixel BC>",
		     80, -0.2, 0.2 );
  TProfile2D * pixqvsdxdy = new
    TProfile2D( "pixqvsdxdy",
		"DUT pixel signal map;height [mm];depth [mm];<pixel signal> [ToT]",
		80, -0.200, 0.200, 60, -0.150, 0.150 );

  TH1I pixdycHisto( "pixdyc",
		    "pixel - Telescope dy;pixel - triplet #Deltaxy [mm];pixels",
		    200, -0.5, 0.5 );
  TProfile pixdyvsz( "pixdyvsz",
		     "#Deltay vs z;z [mm];<pixel - triplet #Deltay> [mm]",
		     200, -10, 10, -0.1, 0.1 );
  TProfile pixdyvsy( "pixdyvsy",
		     "#Deltay vs y;track y [mm];<pixel - triplet #Deltay> [mm]",
		     100, -5, 5, -0.1, 0.1 );
  TProfile pixdyvsty( "pixdyvsty",
		      "#Deltay vs #theta_{y};track y slope [mrad];<pixel - triplet #Deltay> [#mum]",
		      40, -2, 2, -100, 100 );
  TProfile pixdyvsev( "pixdyvsev", "PIX - Telescope y vs time;time [events];<#Deltay> [mm]",
		      600, 0, 600*1000, -0.1, 0.1 );

  TH1I roadnpxHisto( "roadnpx",
		    "pixels in track road;in-road pixels;fiducial tracks",
		    400,  0.5, 400.5 );
  TH1I roadncolHisto( "roadncol",
		      "filled columns in track road;filled columns;tracks",
		      150, 0.5, 150.5 );
  TProfile roadncolvscol0( "roadncolvscol0",
			   "filled columns in track road;first column;filled columns;tracks",
			   400, 0, 400 );

  TH1I tridHisto( "trid",
		  "in-time triplet depth;depth [mm];in-time triplets",
		  80, -0.2, 0.2 );
  TH1I tridlkHisto( "tridlk",
		  "in-time triplet depth;depth [mm];in-time triplet hits",
		  80, -0.2, 0.2 );
  TProfile pixqvsd( "pixqvsd",
		    "PIX column signal vs depth [mm];depth [mm];<column signal> [ToT]",
		    80, -0.2, 0.2 );

  TH1I triymHisto( "triym",
		   "in-time in-depth triplet y;y [mm];in-time in-depth triplets",
		   120, -6, 6 );
  TH1I triymlkHisto( "triymlk",
		     "in-time in-depth linked triplet y;y [mm];in-time in-depth triplet hits",
		     120, -6, 6 );

  TProfile modlkvst1( "modlkvst1",
		      "triplet-MOD links vs time;time [s];triplets with MOD links / s",
		      300, 0, 300, -0.5, 1.5 );
  TProfile modlkvst3( "modlkvst3",
		      "triplet-MOD links vs time;time [s];triplets with MOD links / 10s",
		      150, 0, 1500, -0.5, 1.5 );
  TProfile modlkvst5( "modlkvst5",
		      "triplet-MOD links vs time;time [s];triplets with MOD links / min",
		      1100, 0, 66000, -0.5, 1.5 );
  TProfile modlkvsev( "modlkvsev",
		      "triplet-MOD links vs events;events;triplets with MOD links / 1000",
		      1100, 0, 1.1E6, -0.5, 1.5 );

  TH1I ntrimodHisto( "ntrimod", "triplet - MOD links;triplet - MOD links;events",
		    11, -0.5, 10.5 );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // event loop:

  cout << endl;

  FileReader * reader;
  if(      run <    100 )
    reader = new FileReader( runnum.c_str(), "data/run0000$2R$X" );
  else if( run <   1000 )
    reader = new FileReader( runnum.c_str(), "data/run000$3R$X" );
  else if( run <  10000 )
    reader = new FileReader( runnum.c_str(), "data/run00$4R$X" );
  else if( run < 100000 )
    reader = new FileReader( runnum.c_str(), "data/run0$5R$X" );
  else
    reader = new FileReader( runnum.c_str(), "data/run$6R$X" );

  DetectorEvent evt = reader->GetDetectorEvent();
  uint64_t evTLU0 = evt.GetTimestamp(); // 384 MHz = 2.6 ns
  if( evt.IsBORE() ) {
    eudaq::PluginManager::Initialize(evt);
    cout << "BORE TLU " << evTLU0 << endl;
    reader->NextEvent();
  }
  uint64_t prevTLU = evTLU0;

  int iev = 0;

  if( fev ) cout << "EU skip " << fev << endl;
  while( iev < fev ) {
    reader->NextEvent(); // fast forward
    ++iev;
  }

  int nevA = 0;
  int nevB = 0;

  do {
  
    evt = reader->GetDetectorEvent();

    uint64_t evTLU = evt.GetTimestamp(); // 384 MHz = 2.6 ns

    double evsec = (evTLU - evTLU0) / fTLU;
    t1Histo.Fill( evsec );
    t2Histo.Fill( evsec );
    t3Histo.Fill( evsec );
    t4Histo.Fill( evsec );
    t5Histo.Fill( evsec );
    t6Histo.Fill( evsec/3600 );

    double evdt = (evTLU - prevTLU) / fTLU;
    dtusHisto.Fill( evdt * 1E6 ); // [us]
    dtmsHisto.Fill( evdt * 1E3 ); // [ms]

    dt373Histo.Fill( (evTLU - prevTLU)%373 );
    dt374Histo.Fill( (evTLU - prevTLU)%374 );
    dt375Histo.Fill( (evTLU - prevTLU)%375 ); // best
    dt376Histo.Fill( (evTLU - prevTLU)%376 );
    dt377Histo.Fill( (evTLU - prevTLU)%377 );
    dt375vsdt.Fill( evdt*1E6, (evTLU - prevTLU)%375 ); // linear

    prevTLU = evTLU;

    bool ldbg = 0;

    if( iev <  0 )
      ldbg = 1;

    if( lev < 100 )
      ldbg = 1;

    if( iev < 10 || ldbg )
      cout << "edg53 processing  " << run << "." << iev << "  taken " << evsec << endl;
    else if( iev < 100 && iev%10 == 0 )
      cout << "edg53 processing  " << run << "." << iev << "  taken " << evsec << endl;
    else if( iev < 1000 && iev%100 == 0 )
      cout << "edg53 processing  " << run << "." << iev << "  taken " << evsec << endl;
    else if( iev%1000 == 0 ) {
      cout << "edg53 processing  " << run << "." << iev
	   << "  taken " << evsec
	   << endl;
    }

    StandardEvent sevt = eudaq::PluginManager::ConvertToStandard(evt);

    if( ldbg ) cout << "planes " << sevt.NumPlanes() << endl;

    vector <pixel> pbDUT;
    vector < cluster > cl[9];

    for( size_t iplane = 0; iplane < sevt.NumPlanes(); ++iplane ) {

      const eudaq::StandardPlane &plane = sevt.GetPlane(iplane);

      if( ldbg )
	std::cout
	  << "  " << iplane
	  << ": plane " << plane.ID() // 1
	  << " " << plane.Type() // NI or BDAQ53
	  << " " << plane.Sensor() // MIMOSA26 or RD53A
	  << " frames " << plane.NumFrames() // 2 for NI or 32 for BDAQ53
	  << " pivot " << plane.PivotPixel()
	  << " total " << plane.TotalPixels()
	  << " hits " << plane.HitPixels()
	  ;

      int ipl = plane.ID(); // 0 = DUT, 1..6 = Mimosa

      if( ipl < 0 || ipl > 6 ) {
	cout << "event " << iev << " wrong plane number " << ipl << endl;
	continue;
      }

      hpivot[ipl].Fill( plane.PivotPixel() );
      hnpx[ipl].Fill( plane.HitPixels() );
      hnframes[ipl].Fill( plane.NumFrames() ); // 32
      if( ipl == iDUT )
	dutnpxvsev.Fill( iev, plane.HitPixels() );

      vector <pixel> pb; // for clustering

      // loop over frames, then pixels per frame

      for( unsigned frm = 0; frm < plane.NumFrames(); ++frm ) 

	for( size_t ipix = 0; ipix < plane.HitPixels( frm ); ++ipix ) {

	  if( ldbg ) 
	    std::cout << ": " << plane.GetX(ipix,frm)
		      << "." << plane.GetY(ipix,frm)
		      << "." << plane.GetPixel(ipix,frm) << " ";

	  int ix = plane.GetX(ipix,frm); // column
	  int iy = plane.GetY(ipix,frm); // row
	  int tot = plane.GetPixel(ipix,frm); // ToT 0..15

	  if( ipl == iDUT ) {
	    dutpxbc0Histo.Fill( frm ); // before hot pixel masking
	    hmap[8]->Fill( ix+0.5, iy+0.5 ); // before masking
	  }

	  // skip hot pixels:

	  int ipx = ix*ny[ipl] + iy;
	  if( hotset[ipl].count(ipx) ) continue;

	  pixel px;
	  px.col = ix; // ROC col
	  px.row = iy; // row
	  px.tot = tot;
	  px.frm = frm;

	  if( ipl == iDUT ) {

	    px.tot += 1; // shift from zero

	    if( !fifty ) { // 100x25 from ROC to sensor:
	      px.col = ix/2; // 100 um
	      if( ix%2 ) 
		px.row = 2*iy + 0; // different from R4S
	      else
		px.row = 2*iy + 1; // see ed53 for shallow angle
	      if( chip0 == 182 || chip0 == 211 || chip0 == 512 ) { // HLL
		if( ix%2 ) 
		  px.row = 2*iy + 1;
		else
		  px.row = 2*iy + 0;
	      }
	    }
	  } // DUT

	  pb.push_back(px);

	  hcol[ipl].Fill( ix+0.5 ); // ROC
	  hrow[ipl].Fill( iy+0.5 );
	  hmap[ipl]->Fill( ix+0.5, iy+0.5 ); // after thr

	} // pix

      if( ldbg ) std::cout << std::endl;

      // clustering:

      if( ipl == iDUT )
	pbDUT = pb; // no clustering
      else
	cl[ipl] = getClusn( pb ); // Mimosa

      if( ldbg ) cout << "    clusters " << cl[ipl].size() << endl;

      hncl[ipl].Fill( cl[ipl].size() );

      for( vector<cluster>::iterator c = cl[ipl].begin(); c != cl[ipl].end(); ++c ) {

	hsiz[ipl].Fill( c->size );
	hncol[ipl].Fill( c->ncol );
	hnrow[ipl].Fill( c->nrow );

	// cluster isolation:

	vector<cluster>::iterator d = c; // upper diagonal
	++d;
	for( ; d != cl[ipl].end(); ++d ) {
	  double dx = d->col - c->col;
	  double dy = d->row - c->row;
	  double dxy = sqrt( dx*dx + dy*dy );
	  if( dxy < c->mindxy ) c->mindxy = dxy;
	  if( dxy < d->mindxy ) d->mindxy = dxy;
	}

      } // cl

      for( vector<cluster>::iterator c = cl[ipl].begin(); c != cl[ipl].end(); ++c )
	hdxy[ipl].Fill( c->mindxy );

    } // eudaq planes

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // DUT:

    for( unsigned ipx = 0; ipx < pbDUT.size(); ++ipx ) {

      int col = pbDUT[ipx].col; // sensor
      int row = pbDUT[ipx].row; // sensor
      int tot = pbDUT[ipx].tot;
      int frm = pbDUT[ipx].frm; // [BC]

      dutpxcolHisto.Fill( col + 0.5 );
      dutpxrowHisto.Fill( row + 0.5 );
      dutpxqHisto.Fill( tot + 0.5 );
      dutpxbcHisto.Fill( frm + 0.5 );

      dutpxqvsx.Fill( col + 0.5, tot + 0.5 );
      dutpxqvsxy->Fill( col + 0.5, row + 0.5, tot + 0.5 );

    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // MOD:

    if( iev >= fev && modrun &&
	Astream.good() && ! Astream.eof() &&
	Bstream.good() && ! Bstream.eof()
	) {

      // one line = one trigger from one TBM channel

      string sl;
      getline( Astream, sl );

      istringstream Aev( sl ); // tokenize string

      //687
      //688 265 74 79 266 73 74 266 74 66

      int trg;
      Aev >> trg;

      ++nevA;

      int ierr = 0;

      bool ldb = 0;

      vector <pixel> pbmod; // for clustering

      while( ! Aev.eof()
	     // && Aev.good() && Aev.tellg() > 0 && Aev.tellg() < (int) sl.size()
	     ) {

	int xm;
	Aev >> xm; // 0..415

	if( Aev.eof() ) {
	  cout << "A truncated at line " << nevA << " col " << xm << endl;
	  cout << "Aev " << sl << endl;
	  ierr = 1;
	  break;
	}
	if( ldb ) cout << " : " << xm << flush;

	int ym;
	Aev >> ym; // 0..159

	if( Aev.eof() ) {
	  cout << "A truncated at line " << nevA << " row " << ym << endl;
	  cout << "Aev " << sl << endl;
	  ierr = 1;
	  break;
	}
	if( ldb ) cout << "." << ym << flush;

	int adc;
	Aev >> adc;
	if( ldb ) cout << "." << adc << flush;
	if( ldb ) cout << " (" << Aev.tellg() << ") " << flush;

	if( xm < 0 || xm > 415 || ym < 0 || ym > 159 || adc < 0 || adc > 255 ) {
	  cout << "data error at line " << nevA << " ev " << trg << endl;
	  ierr = 1;
	  break;
	}

	int roc = xm / 52; // 0..7
	int col = xm % 52; // 0..51
	//int row = ym;

	// leave space for big pixels:

	int ix = 1 + xm + 2*roc; // 1..52 per ROC
	int iy = ym;
	if( ym > 79 ) iy += 2; // 0..79, 82..161

	// fill pixel block for clustering:

	pixel px;
	px.col = ix; // col
	px.row = iy; // row
	px.tot = adc;
	pbmod.push_back(px);

	// double big pixels:
	// 0+1
	// 2..51
	// 52+53

	if( col == 0 ) {
	  px.col = ix-1; // double
	  px.row = iy;
	  pbmod[pbmod.size()-1].tot *= 0.5;
	  px.tot = 0.5*adc;
	  pbmod.push_back(px);
	}

	if( col == 51 ) {
	  px.col = ix+1; // double
	  px.row = iy;
	  pbmod[pbmod.size()-1].tot *= 0.5;
	  px.tot = 0.5*adc;
	  pbmod.push_back(px);
	}

	if( ym == 79 ) {
	  px.col = ix;
	  px.row = 80; // double
	  pbmod[pbmod.size()-1].tot *= 0.5;
	  px.tot = 0.5*adc;
	  pbmod.push_back(px);
	}

	if( ym == 80 ) {
	  px.col = ix;
	  px.row = 81; // double
	  pbmod[pbmod.size()-1].tot *= 0.5;
	  px.tot = 0.5*adc;
	  pbmod.push_back(px);
	}

      } // Aev

      if( ldb ) cout << endl;

      // B:

      getline( Bstream, sl );

      istringstream Bev( sl ); // tokenize string

      // 97875 358  81 47
      // 98144 373 126 47 372 126 46

      Bev >> trg;

      ++nevB; // lines

      if( ldb ) cout << trg << " B";

      //if( Bev.eof() ) { cout << "B empty" << endl; continue; }

      while( ! Bev.eof()
	     //&& Bev.good() && Bev.tellg() > 0 && Bev.tellg() < (int) sl.size()
	     ) {

	int xm;
	Bev >> xm; // 0..415

	if( Bev.eof() ) {
	  cout << "B truncated at line " << nevB << " col " << xm << endl;
	  cout << "Bev " << sl << endl;
	  ierr = 1;
	  break;
	}
	if( ldb ) cout << " : " << xm << flush;

	int ym;
	Bev >> ym; // 0..159

	if( Bev.eof() ) {
	  cout << "B truncated at line " << nevB << " row " << ym << endl;
	  cout << "Bev " << sl << endl;
	  ierr = 1;
	  break;
	}
	if( ldb ) cout << "." << ym << flush;

	int adc;
	Bev >> adc;
	if( ldb ) cout << "." << adc << flush;
	if( ldb ) cout << " (" << Bev.tellg() << ") " << flush;

	if( xm < 0 || xm > 415 || ym < 0 || ym > 159 || adc < 0 || adc > 255 ) {
	  cout << "data error at line " << nevB << " ev " << trg << endl;
	  ierr = 1;
	  break;
	}

	int roc = xm / 52; // 0..7
	int col = xm % 52; // 0..51
	//int row = ym;

	// leave space for big pixels:

	int ix = 1 + xm + 2*roc; // 1..52 per ROC
	int iy = ym;
	if( ym > 79 ) iy += 2; // 0..79, 82..161

	// fill pixel block for clustering:

	pixel px;
	px.col = ix; // col
	px.row = iy; // row
	px.tot = adc;
	pbmod.push_back(px);

	// double big pixels:
	// 0+1
	// 2..51
	// 52+53

	if( col == 0 ) {
	  px.col = ix-1; // double
	  px.row = iy;
	  pbmod[pbmod.size()-1].tot *= 0.5;
	  px.tot = 0.5*adc;
	  pbmod.push_back(px);
	}

	if( col == 51 ) {
	  px.col = ix+1; // double
	  px.row = iy;
	  pbmod[pbmod.size()-1].tot *= 0.5;
	  px.tot = 0.5*adc;
	  pbmod.push_back(px);
	}

	if( ym == 79 ) {
	  px.col = ix;
	  px.row = 80; // double
	  pbmod[pbmod.size()-1].tot *= 0.5;
	  px.tot = 0.5*adc;
	  pbmod.push_back(px);
	}

	if( ym == 80 ) {
	  px.col = ix;
	  px.row = 81; // double
	  pbmod[pbmod.size()-1].tot *= 0.5;
	  px.tot = 0.5*adc;
	  pbmod.push_back(px);
	}

      } // Bev

      if( ldb ) cout << endl;

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // module clustering:

      hnpx[iMOD].Fill( pbmod.size() );

      if( !ierr )
	cl[iMOD] = getClusq( pbmod );

      hncl[iMOD].Fill( cl[iMOD].size() );

      for( vector<cluster>::iterator c = cl[iMOD].begin(); c != cl[iMOD].end(); ++c ) {

	hsiz[iMOD].Fill( c->size );
	hncol[iMOD].Fill( c->ncol );
	hnrow[iMOD].Fill( c->nrow );

	for( vector<pixel>::iterator px = c->vpix.begin(); px != c->vpix.end(); ++px ) {
	  hcol[iMOD].Fill( px->col+0.5 );
	  hrow[iMOD].Fill( px->row+0.5 );
	  hmap[iMOD]->Fill( px->col+0.5, px->row+0.5 );
	}

      } // cl

    } // mod

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // make triplets 1+3-2:

    vector <triplet> triplets;

    //double triCut = 0.1; // [mm]
    double triCut = 0.05; // [mm] like tele

    for( vector<cluster>::iterator cA = cl[1].begin(); cA != cl[1].end(); ++cA ) {

      double xA = cA->col*ptchx[1] - alignx[1];
      double yA = cA->row*ptchy[1] - aligny[1];
      double xmid = xA - midx[1];
      double ymid = yA - midy[1];
      xA = xmid - ymid*rotx[1];
      yA = ymid + xmid*roty[1];

      double zA = zz[1] + alignz[1];
      double zC = zz[3] + alignz[3];
      double zB = zz[2] + alignz[2];

      for( vector<cluster>::iterator cC = cl[3].begin(); cC != cl[3].end(); ++cC ) {

	double xC = cC->col*ptchx[3] - alignx[3];
	double yC = cC->row*ptchy[3] - aligny[3];
	double xmid = xC - midx[3];
	double ymid = yC - midy[3];
	xC = xmid - ymid*rotx[3];
	yC = ymid + xmid*roty[3];

	double dx2 = xC - xA;
	double dy2 = yC - yA;
	double dzCA = zC - zA;
	hdx13.Fill( dx2 );
	hdy13.Fill( dy2 );

	if( fabs( dx2 ) > 0.005*f * dzCA ) continue; // angle cut *f?
	if( fabs( dy2 ) > 0.005*f * dzCA ) continue; // angle cut

	double avx = 0.5 * ( xA + xC ); // mid
	double avy = 0.5 * ( yA + yC );
	double avz = 0.5 * ( zA + zC ); // mid z
 
	double slpx = ( xC - xA ) / dzCA; // slope x
	double slpy = ( yC - yA ) / dzCA; // slope y

	// middle plane B = 2:

	for( vector<cluster>::iterator cB = cl[2].begin(); cB != cl[2].end(); ++cB ) {

	  double xB = cB->col*ptchx[2] - alignx[2];
	  double yB = cB->row*ptchy[2] - aligny[2];
	  double xmid = xB - midx[2];
	  double ymid = yB - midy[2];
	  xB = xmid - ymid*rotx[2];
	  yB = ymid + xmid*roty[2];

	  // interpolate track to B:

	  double dz = zB - avz;
	  double xm = avx + slpx * dz; // triplet at mid
	  double ym = avy + slpy * dz;

	  double dxm = xB - xm;
	  double dym = yB - ym;
	  htridx.Fill( dxm );
	  htridy.Fill( dym );

	  if( fabs( dym ) < 0.05 ) {

	    htridxc.Fill( dxm );
	    tridxvsx.Fill( xB, dxm );
	    tridxvsy.Fill( yB, dxm );
	    tridxvstx.Fill( slpx, dxm );
	    tridxvst3.Fill( evsec, dxm );
	    tridxvst5.Fill( evsec, dxm );

	  } // dy

	  if( fabs( dxm ) < 0.05 ) {
	    htridyc.Fill( dym );
	    tridyvsx.Fill( xB, dym );
	    tridyvsty.Fill( slpy, dym );
	    tridyvst3.Fill( evsec, dym );
	    tridyvst5.Fill( evsec, dym );
	  }

	  // telescope triplet cuts:

	  if( fabs(dxm) > triCut ) continue;
	  if( fabs(dym) > triCut ) continue;

	  triplet tri;
	  tri.xm = avx;
	  tri.ym = avy;
	  tri.zm = avz;
	  tri.sx = slpx;
	  tri.sy = slpy;
	  tri.lk = 0;
	  tri.ttdmin = 99.9; // isolation [mm]
	  tri.iA = distance( cl[1].begin(), cA );
	  tri.iB = distance( cl[2].begin(), cB );
	  tri.iC = distance( cl[3].begin(), cC );

	  triplets.push_back(tri);

	  trixHisto.Fill( avx );
	  triyHisto.Fill( avy );
	  trixyHisto->Fill( avx, avy );
	  tritxHisto.Fill( slpx );
	  trityHisto.Fill( slpy );

	} // cl B

      } // cl C

    } // cl A

    ntriHisto.Fill( triplets.size() );

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // triplets vs MOD and DUT:

    int nmdm = 0;
    int ntrimod = 0;

    double xcutMOD = 0.15;
    double ycutMOD = 0.15;

    for( unsigned int iA = 0; iA < triplets.size(); ++iA ) { // iA = upstream

      double xmA = triplets[iA].xm;
      double ymA = triplets[iA].ym;
      double zmA = triplets[iA].zm;
      double sxA = triplets[iA].sx;
      double syA = triplets[iA].sy;

      double zB = MODz - zmA; // z MOD from mid of triplet
      double xB = xmA + sxA * zB; // triplet impact point on MOD
      double yB = ymA + syA * zB;

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // tri vs tri: isolation at MOD

      double ttdminmod = 99.9;

      for( unsigned int jj = 0; jj < triplets.size(); ++jj ) {

	if( jj == iA ) continue;

	double xmj = triplets[jj].xm;
	double ymj = triplets[jj].ym;
	double sxj = triplets[jj].sx;
	double syj = triplets[jj].sy;

	double dz = MODz - triplets[jj].zm;
	double xj = xmj + sxj * dz; // triplet impact point on MOD
	double yj = ymj + syj * dz;

	double dx = xB - xj;
	double dy = yB - yj;
	double dd = sqrt( dx*dx + dy*dy );
	if( dd < ttdminmod )
	  ttdminmod = dd;

      } // jj

      ttdminmod1Histo.Fill( ttdminmod );
      ttdminmod2Histo.Fill( ttdminmod );
      triplets[iA].ttdmin = ttdminmod;

      // intersect inclined track with tilted MOD plane:

      double zd = (Nzm*zB - Nym*ymA - Nxm*xmA) / (Nxm*sxA + Nym*syA + Nzm); // from zmB
      double yd = ymA + syA * zd;
      double xd = xmA + sxA * zd;

      double dzd = zd + zmA - MODz; // from MOD z0 [-8,8] mm

      // transform into MOD system: (passive).
      // large rotations don't commute: careful with order

      double x1m = com*xd - som*dzd; // turn o
      double y1m = yd;
      double z1m = som*xd + com*dzd;

      double x2m = x1m;
      double y2m = cam*y1m + sam*z1m; // tilt a

      double x3m = cfm*x2m + sfm*y2m; // rot
      double y3m =-sfm*x2m + cfm*y2m;

      double x4m =-x3m + MODalignx; // shift to mid
      double y4m = y3m + MODaligny; // invert y, shift to mid

      double xmodm = fmod( 36.000 + x4m, 0.3 ); // [0,0.3] mm, 2 pixel wide
      double ymodm = fmod(  9.000 + y4m, 0.2 ); // [0,0.2] mm

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // triplets vs MOD clusters:

      bool ltrimod = 0;

      for( vector<cluster>::iterator c = cl[iMOD].begin(); c != cl[iMOD].end(); ++c ) {

	double ccol = c->col;
	double crow = c->row;
	double modx = ( ccol + 0.5 - nx[iMOD]/2 ) * ptchx[iMOD]; // -33..33 mm
	double mody = ( crow + 0.5 - ny[iMOD]/2 ) * ptchy[iMOD]; // -8..8 mm
	double q = c->signal;
	double q0 = q*normm;

	modxHisto.Fill( modx );
	modyHisto.Fill( mody );

	// residuals for pre-alignment:

	modsxaHisto.Fill( modx + x3m ); // peak
	moddxaHisto.Fill( modx - x3m ); // 

	modsyaHisto.Fill( mody + y3m ); // 
	moddyaHisto.Fill( mody - y3m ); // peak

	double moddx = modx - x4m;
	double moddy = mody - y4m;

	moddxHisto.Fill( moddx );
	moddyHisto.Fill( moddy );

	if( fabs( moddx ) < xcutMOD ) {

	  moddycHisto.Fill( moddy );
	  moddyvsx.Fill( -x3m, moddy ); // for rot
	  moddyvsy.Fill( y2m, moddy ); // for tilt
	  moddyvsty.Fill( syA, moddy );
	  moddyvst5.Fill( evsec, moddy );

	}

	if( fabs( moddy ) < ycutMOD ) {

	  moddxcHisto.Fill( moddx );
	  moddxvsx.Fill( -x1m, moddx ); // for turn
	  moddxvsy.Fill( y3m, moddx ); // for rot
	  moddxvstx.Fill( sxA, moddx );
	  moddxvst5.Fill( evsec, moddx );

	}

	if( fabs( moddx ) < xcutMOD &&
	    fabs( moddy ) < ycutMOD ) {

	  modnpxHisto.Fill( c->size );
	  modqHisto.Fill( q );
	  modq0Histo.Fill( q0 );
	  modnpxvsxmym->Fill( xmodm*1E3, ymodm*1E3, c->size );

	}

	if( fabs( moddx ) < xcutMOD &&
	    fabs( moddy ) < ycutMOD ) {

	  modlkxBHisto.Fill( xB );
	  modlkyBHisto.Fill( yB );
	  modlkxHisto.Fill( x4m );
	  modlkyHisto.Fill( y4m );
	  modlkcolHisto.Fill( ccol );
	  modlkrowHisto.Fill( crow );

	  //if( crow > 80 ) cout << "B link " << crow << " ev " << iev << endl;
	  //if( crow < 80 ) cout << "A link " << crow << " ev " << iev << endl;
	  //if( iev > 5100200 && crow < 80 ) cout << "link " << crow << " ev " << iev << endl; // A

	  triplets[iA].lk = 1;
	  nmdm = 1; // we have a MOD-triplet match in this event
	  ++ntrimod;
	  ltrimod = 1;

	} // MOD link x and y

      } // MOD

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // tri vs tri: isolation at DUT

      double dzA = DUTz - zmA; // z DUT from mid of triplet
      double xAc = xmA + sxA * dzA; // track at z_mid(DUT)
      double yAc = ymA + syA * dzA;

      double ttdmin = 99.9;

      for( unsigned int jj = 0; jj < triplets.size(); ++jj ) {

	if( jj == iA ) continue;

	double xmj = triplets[jj].xm;
	double ymj = triplets[jj].ym;
	double sxj = triplets[jj].sx;
	double syj = triplets[jj].sy;

	double dz = dzA + zmA - triplets[jj].zm;
	double xj = xmj + sxj * dz; // triplet impact point on DUT
	double yj = ymj + syj * dz;

	double dx = xAc - xj;
	double dy = yAc - yj;
	double dd = sqrt( dx*dx + dy*dy );
	if( dd < ttdmin )
	  ttdmin = dd;

      } // jj

      ttdmin1Histo.Fill( ttdmin );
      ttdmin2Histo.Fill( ttdmin );
      triplets[iA].ttdmin = ttdmin;

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // triplet vs DUT:

      if( ltrimod ) { // in-time

	trixcHisto.Fill( xAc );
	triycHisto.Fill( yAc );

	int nroad = 0;
	vector <int> roadcol(nx[iDUT]);

	vector <int> colq(nx[iDUT]);

	for( unsigned ipx = 0; ipx < pbDUT.size(); ++ipx ) { // pixels

	  int col = pbDUT[ipx].col; // sensor
	  int row = pbDUT[ipx].row; // sensor
	  int tot = pbDUT[ipx].tot;
	  int frm = pbDUT[ipx].frm; // [BC]

	  double pz = ( col + 0.5 - nx[iDUT]/2 ) * ptchx[iDUT]; // mm
	  double py = ( row + 0.5 - ny[iDUT]/2 ) * ptchy[iDUT]; // mm

	  pixzyHisto->Fill( pz, py );

	  // transform pixel into telescope system:

	  double x1 = 0; // sensor mid plane
	  double y1 = cf*py + sf*pz; // rot
	  double z1 =-sf*py + cf*pz;

	  double x2 = co*x1 - so*z1; // turn x-z
	  double y2 = y1;
	  double z2 = so*x1 + co*z1;

	  double x3 = ca*x2 + sa*y1; // tilt y-x
	  double y3 =-sa*x2 + ca*y2;
	  double z3 = z2;

	  // track at z3:

	  double dz = z3 + DUTz - zmA; // z of pixel from mid of triplet A
	  double xA = xmA + sxA * dz; // triplet impact point on DUT
	  double yA = ymA + syA * dz; // track A at pixel

	  pixdxaHisto.Fill( xA - x3 );
	  pixsxaHisto.Fill( xA + x3 );
	  pixdyaHisto.Fill( yA - y3 );
	  pixsyaHisto.Fill( yA + y3 );

	  double dx = xA - x3 - DUTalignx;
	  double dy = yA - y3 - DUTaligny;

	  pixdyHisto.Fill( dy );
	  pixqvsdxdy->Fill( dy, dx, tot );

	  if( fabs( dy ) < 0.07 ) { // track road

	    pixdxcHisto.Fill( dx );
	    pixdxvsz.Fill( z3, dx ); // tan(turn) = slope
	    pixdxvsy.Fill( yA, dx ); // tan(turn) = slope
	    pixdxvsev.Fill( iev, dx );

	    colq[col] += tot;

	    if( fabs( dx ) < 0.150 ) { // depth

	      hroadmap->Fill( col+0.5, row+0.5 );
	      ++nroad;
	      ++roadcol[col];
	      pixbclkHisto.Fill( frm ); // linked
	      pixbcvsd.Fill( dx, frm ); // timewalk? drift?

	    }

	  } // dy

	  if( fabs( dx ) < 0.150 ) { // depth

	    pixdycHisto.Fill( dy ); // aligny
	    pixdyvsz.Fill( z3, dy ); // rot
	    pixdyvsy.Fill( yA, dy ); // tan(tilt) = slope
	    pixdyvsty.Fill( syA*1E3, dy*1E3 ); // slope = -dz
	    pixdyvsev.Fill( iev, dy );

	  }

	} // loop pix

	roadnpxHisto.Fill( nroad );

	int ncol = 0;
	int col0 = nx[iDUT];
	int col9 = 0;
	for( int icol = 0; icol < nx[iDUT]; ++icol ) {
	  if( roadcol[icol] ) {
	    ++ncol;
	    if( icol < col0 ) col0 = icol;
	    if( icol > col9 ) col9 = icol;
	  }
	}

	roadncolHisto.Fill( ncol );
	roadncolvscol0.Fill( col0+0.5, ncol ); // overflows have weight zero

	// get mean depth per column:

	for( int col = 128/(2-fifty); col < 264/(2-fifty); ++col ) { // Lin

	  double depth = 0;
	  int nd = 0;
	  bool in = 0;
	  double ym = yAc; // track y at DUTz
	  double dymin = 999;

	  for( int row = 0; row < 384/(2-fifty); ++row ) {

	    double pz = ( col + 0.5 - nx[iDUT]/2 ) * ptchx[iDUT]; // mm
	    double py = ( row + 0.5 - ny[iDUT]/2 ) * ptchy[iDUT]; // mm

	    // transform pixel into telescope system:

	    double x1 = 0; // sensor mid plane
	    double y1 = cf*py + sf*pz; // rot
	    double z1 =-sf*py + cf*pz;

	    double x2 = co*x1 - so*z1; // turn x-z
	    double y2 = y1;
	    double z2 = so*x1 + co*z1;

	    double x3 = ca*x2 + sa*y1; // tilt y-x
	    double y3 =-sa*x2 + ca*y2;
	    double z3 = z2;

	    // track at z3:

	    double dz = z3 + DUTz - zmA; // z of pixel from mid of triplet A
	    double xA = xmA + sxA * dz; // triplet impact point on DUT
	    double yA = ymA + syA * dz; // track A at pixel

	    double dx = xA - x3 - DUTalignx;
	    double dy = yA - y3 - DUTaligny;

	    if( fabs( dx ) > 0.220 ) continue; // depth
	    if( fabs( dy ) > 0.100 ) continue; // road

	    depth += dx;
	    ++nd;

	    if( fabs( dx ) < 0.075 )
	      in = 1;

	    if( fabs(dy) < fabs(dymin) ) {
	      dymin = dy;
	      ym = yA; // track at pixel
	    }

	  } // rows

	  if( nd ) { // track in depth

	    depth /= nd;

	    tridHisto.Fill( depth );
	    if( roadcol[col] )
	      tridlkHisto.Fill( depth );

	    pixqvsd.Fill( depth, colq[col] );

	  }

	  if( in ) { // track in depth: fiducial

	    triymHisto.Fill( ym );
	    if( roadcol[col] )
	      triymlkHisto.Fill( ym ); // check y acceptance

	  }

	} // cols = z along track

      } // in-time

    } // loop triplets

    modlkvst1.Fill( evsec, nmdm ); // MOD yield vs time
    modlkvst3.Fill( evsec, nmdm );
    modlkvst5.Fill( evsec, nmdm );
    modlkvsev.Fill( iev, nmdm );
    ntrimodHisto.Fill( ntrimod );

    ++iev;

  } while( reader->NextEvent() && iev < lev );

  delete reader;

  cout << "done after " << iev << " events" << endl;
  histoFile->Write();
  histoFile->Close();

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // MOD alignment:

  if( moddxaHisto.GetEntries() > 9999 ) {

    double newMODalignx = MODalignx;
    double newMODaligny = MODaligny;

    if( moddxaHisto.GetMaximum() > modsxaHisto.GetMaximum() ) {
      cout << endl << moddxaHisto.GetTitle()
	   << " bin " << moddxaHisto.GetBinWidth(1)
	   << endl;
      TF1 * fgp0 = new TF1( "fgp0", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -1, 1 );
      double xpk = moddxaHisto.GetBinCenter( moddxaHisto.GetMaximumBin() );
      fgp0->SetParameter( 0, moddxaHisto.GetMaximum() ); // amplitude
      fgp0->SetParameter( 1, xpk );
      fgp0->SetParameter( 2, moddxaHisto.GetBinWidth(1) ); // sigma
      fgp0->SetParameter( 3, moddxaHisto.GetBinContent( moddxaHisto.FindBin(xpk-1) ) ); // BG
      moddxaHisto.Fit( "fgp0", "q", "", xpk-1, xpk+1 );
      cout << "Fit Gauss + BG:"
	   << endl << "  A " << fgp0->GetParameter(0)
	   << endl << "mid " << fgp0->GetParameter(1)
	   << endl << "sig " << fgp0->GetParameter(2)
	   << endl << " BG " << fgp0->GetParameter(3)
	   << endl;
      newMODalignx = fgp0->GetParameter(1);
      delete fgp0;
    }
    else {
      cout << endl << modsxaHisto.GetTitle()
	   << " bin " << modsxaHisto.GetBinWidth(1)
	   << endl;
      TF1 * fgp0 = new TF1( "fgp0", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -1, 1 );
      double xpk = modsxaHisto.GetBinCenter( modsxaHisto.GetMaximumBin() );
      fgp0->SetParameter( 0, modsxaHisto.GetMaximum() ); // amplitude
      fgp0->SetParameter( 1, xpk );
      fgp0->SetParameter( 2, modsxaHisto.GetBinWidth(1) ); // sigma
      fgp0->SetParameter( 3, modsxaHisto.GetBinContent( modsxaHisto.FindBin(xpk-1) ) ); // BG
      modsxaHisto.Fit( "fgp0", "q", "", xpk-1, xpk+1  );
      cout << "Fit Gauss + BG:"
	   << endl << "  A " << fgp0->GetParameter(0)
	   << endl << "mid " << fgp0->GetParameter(1)
	   << endl << "sig " << fgp0->GetParameter(2)
	   << endl << " BG " << fgp0->GetParameter(3)
	   << endl;
      newMODalignx = fgp0->GetParameter(1);
      delete fgp0;
    }

    if( moddyaHisto.GetMaximum() > modsyaHisto.GetMaximum() ) {
      cout << endl << moddyaHisto.GetTitle()
	   << " bin " << moddyaHisto.GetBinWidth(1)
	   << endl;
      TF1 * fgp0 = new TF1( "fgp0", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -1, 1 );
      double xpk = moddyaHisto.GetBinCenter( moddyaHisto.GetMaximumBin() );
      fgp0->SetParameter( 0, moddyaHisto.GetMaximum() ); // amplitude
      fgp0->SetParameter( 1, xpk );
      fgp0->SetParameter( 2, moddyaHisto.GetBinWidth(1) ); // sigma
      fgp0->SetParameter( 3, moddyaHisto.GetBinContent( moddyaHisto.FindBin(xpk-1) ) ); // BG
      moddyaHisto.Fit( "fgp0", "q", "", xpk-1, xpk+1 );
      cout << "Fit Gauss + BG:"
	   << endl << "  A " << fgp0->GetParameter(0)
	   << endl << "mid " << fgp0->GetParameter(1)
	   << endl << "sig " << fgp0->GetParameter(2)
	   << endl << " BG " << fgp0->GetParameter(3)
	   << endl;
      newMODaligny = fgp0->GetParameter(1);
      delete fgp0;
    }
    else {
      cout << endl << modsyaHisto.GetTitle()
	   << " bin " << modsyaHisto.GetBinWidth(1)
	   << endl;
      TF1 * fgp0 = new TF1( "fgp0", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -1, 1 );
      double xpk = modsyaHisto.GetBinCenter( modsyaHisto.GetMaximumBin() );
      fgp0->SetParameter( 0, modsyaHisto.GetMaximum() ); // amplitude
      fgp0->SetParameter( 1, xpk );
      fgp0->SetParameter( 2, modsyaHisto.GetBinWidth(1) ); // sigma
      fgp0->SetParameter( 3, modsyaHisto.GetBinContent( modsyaHisto.FindBin(xpk-1) ) ); // BG
      modsyaHisto.Fit( "fgp0", "q", "", xpk-1, xpk+1 );
      cout << "Fit Gauss + BG:"
	   << endl << "  A " << fgp0->GetParameter(0)
	   << endl << "mid " << fgp0->GetParameter(1)
	   << endl << "sig " << fgp0->GetParameter(2)
	   << endl << " BG " << fgp0->GetParameter(3)
	   << endl;
      newMODaligny = fgp0->GetParameter(1);
      delete fgp0;
    }

    cout << endl << "coarse MODalign x changed by " << newMODalignx - MODalignx << " mm" << endl;

    cout << endl << "coarse MODalign y changed by " << newMODaligny - MODaligny << " mm" << endl;

    // finer alignment:

    if( MODaligniteration > 0 && fabs( newMODalignx - MODalignx ) < 0.1 ) {

      cout << endl << moddxcHisto.GetTitle()
	   << " bin " << moddxcHisto.GetBinWidth(1)
	   << endl;
      TF1 * fgp0 = new TF1( "fgp0", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -1, 1 );
      fgp0->SetParameter( 0, moddxcHisto.GetMaximum() ); // amplitude
      fgp0->SetParameter( 1, moddxcHisto.GetBinCenter( moddxcHisto.GetMaximumBin() ) );
      fgp0->SetParameter( 2, 8*moddxcHisto.GetBinWidth(1) ); // sigma
      fgp0->SetParameter( 3, moddxcHisto.GetBinContent(1) ); // BG
      moddxcHisto.Fit( "fgp0", "q" );
      cout << "Fit Gauss + BG:"
	   << endl << "  A " << fgp0->GetParameter(0)
	   << endl << "mid " << fgp0->GetParameter(1)
	   << endl << "sig " << fgp0->GetParameter(2)
	   << endl << " BG " << fgp0->GetParameter(3)
	   << endl;
      newMODalignx = MODalignx + fgp0->GetParameter(1);
      delete fgp0;

      // dxvsx -> turn:

      if( fabs(som) > 0.01 &&
	  moddxvsx.GetEntries() > 999
	  ) {

	double x0 = -midx[iMOD]+0.2; // fit range
	for( int ix = 1; ix < moddxvsx.GetNbinsX(); ++ix ) {
	  if( moddxvsx.GetBinEntries( ix ) > 11 ) {
	    x0 = moddxvsx.GetBinLowEdge(ix) + 2*moddxvsx.GetBinWidth(ix);
	    break;
	  }
	}

	double x9 = midx[iMOD]-0.2; // [mm] full range
	for( int ix = moddxvsx.GetNbinsX(); ix > 0; --ix ) {
	  if( moddxvsx.GetBinEntries( ix ) > 11 ) {
	    x9 = moddxvsx.GetBinLowEdge(ix)-moddxvsx.GetBinWidth(ix);
	    break;
	  }
	}

	moddxvsx.Fit( "pol1", "q", "", x0, x9 );

	TF1 * fdxvsx = moddxvsx.GetFunction( "pol1" );
	cout << endl << moddxvsx.GetTitle()
	     << ": slope " << fdxvsx->GetParameter(1)
	     << ", extra turn " << fdxvsx->GetParameter(1)/wt/som
	     << " deg"
	     << endl;
	MODturn += fdxvsx->GetParameter(1)/wt/som; // [deg] min 0.6 deg
	//delete fdxvsx;
      }

    } // finer x

    if( MODaligniteration > 0 && fabs( newMODaligny - MODaligny ) < 0.1 ) {

      cout << endl << moddycHisto.GetTitle()
	   << " bin " << moddycHisto.GetBinWidth(1)
	   << endl;
      TF1 * fgp0 = new TF1( "fgp0", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -1, 1 );
      fgp0->SetParameter( 0, moddycHisto.GetMaximum() ); // amplitude
      fgp0->SetParameter( 1, moddycHisto.GetBinCenter( moddycHisto.GetMaximumBin() ) );
      fgp0->SetParameter( 2, 5*moddycHisto.GetBinWidth(1) ); // sigma
      fgp0->SetParameter( 3, moddycHisto.GetBinContent(1) ); // BG
      moddycHisto.Fit( "fgp0", "q" );
      cout << "Fit Gauss + BG:"
	   << endl << "  A " << fgp0->GetParameter(0)
	   << endl << "mid " << fgp0->GetParameter(1)
	   << endl << "sig " << fgp0->GetParameter(2)
	   << endl << " BG " << fgp0->GetParameter(3)
	   << endl;
      newMODaligny = MODaligny + fgp0->GetParameter(1);
      delete fgp0;

      // dyvsx -> rot

      if( moddyvsx.GetEntries() > 999 ) {

	double x0 = -midx[iMOD]+0.2; // fit range
	for( int ix = 1; ix < moddyvsx.GetNbinsX(); ++ix ) {
	  if( moddyvsx.GetBinEntries( ix ) > 11 ) {
	    x0 = moddyvsx.GetBinLowEdge(ix);
	    break;
	  }
	}

	double x9 = midx[iMOD]-0.2;
	for( int ix = moddyvsx.GetNbinsX(); ix > 0; --ix ){
	  if( moddyvsx.GetBinEntries( ix ) > 11 ) {
	    x9 = moddyvsx.GetBinLowEdge(ix)+moddyvsx.GetBinWidth(ix);
	    break;
	  }
	}

	moddyvsx.Fit( "pol1", "q", "", x0, x9 );

	TF1 * fdyvsx = moddyvsx.GetFunction( "pol1" );
	cout << endl << moddyvsx.GetTitle()
	     << ": extra rot " << fdyvsx->GetParameter(1)*1E3 << " mrad" << endl;
	MODrot += fdyvsx->GetParameter(1);
	//delete fdyvsx;

      }

      // dyvsy -> tilt:

      if( fabs( sam ) > 0.01 &&
	  moddyvsy.GetEntries() > 999
	  ) {

	double x0 = -midy[iMOD]+0.2; // fit range
	for( int ix = 1; ix < moddyvsy.GetNbinsX(); ++ix ){
	  if( moddyvsy.GetBinEntries( ix ) > 11 ) {
	    x0 = moddyvsy.GetBinLowEdge(ix);
	    break;
	  }
	}

	double x9 = midy[iMOD]-0.2;
	for( int ix = moddyvsy.GetNbinsX(); ix > 0; --ix ){
	  if( moddyvsy.GetBinEntries( ix ) > 11 ) {
	    x9 = moddyvsy.GetBinLowEdge(ix)+moddyvsy.GetBinWidth(ix);
	    break;
	  }
	}

	moddyvsy.Fit( "pol1", "q", "", x0, x9 );

	TF1 * fdyvsy = moddyvsy.GetFunction( "pol1" );
	cout << endl << moddyvsy.GetTitle()
	     << ": slope " << fdyvsy->GetParameter(1)
	     << ", extra tilt " << fdyvsy->GetParameter(1)/wt/sam
	     << " deg"
	     << endl;

	MODtilt += fdyvsy->GetParameter(1)/wt/sam; // [deg] min 0.6 deg
	//delete fdyvsy;
      }

      // dyvsty -> dz:

      if( moddyvsty.GetEntries() > 999 ) {

	double x0 = -0.002;
	for( int ix = 1; ix < moddyvsty.GetNbinsX(); ++ix ){
	  if( moddyvsty.GetBinEntries( ix ) > 11 ) {
	    x0 = moddyvsty.GetBinLowEdge(ix);
	    break;
	  }
	}

	double x9 = 0.002;
	for( int ix = moddyvsty.GetNbinsX(); ix > 0; --ix ){
	  if( moddyvsty.GetBinEntries( ix ) > 11 ) {
	    x9 = moddyvsty.GetBinLowEdge(ix)+moddyvsty.GetBinWidth(ix);
	    break;
	  }
	}

	moddyvsty.Fit( "pol1", "q", "", x0, x9 );

	TF1 * fdyvsty = moddyvsty.GetFunction( "pol1" );
	cout << endl << moddyvsty.GetTitle()
	     << ": z shift " << fdyvsty->GetParameter(1)
	     << " mm"
	     << endl;
	MODz += fdyvsty->GetParameter(1);
	//delete fdyvsty;
      }

    } // finer y

    // write new MOD alignment:

    ofstream MODalignFile( MODalignFileName.str() );

    MODalignFile << "# MOD alignment for run " << run << endl;
    ++MODaligniteration;
    MODalignFile << "iteration " << MODaligniteration << endl;
    MODalignFile << "alignx " << newMODalignx << endl;
    MODalignFile << "aligny " << newMODaligny << endl;
    MODalignFile << "rot " << MODrot << endl;
    MODalignFile << "tilt " << MODtilt << endl;
    MODalignFile << "turn " << MODturn << endl;
    MODalignFile << "dz " << MODz - zz[1] << endl;

    MODalignFile.close();

    cout << endl << "wrote MOD alignment iteration " << MODaligniteration
	 << " to " << MODalignFileName.str() << endl
	 << "  alignx " << newMODalignx << endl
	 << "  aligny " << newMODaligny << endl
	 << "  rot    " << MODrot << endl
	 << "  tilt   " << MODtilt << endl
	 << "  turn   " << MODturn << endl
	 << "  dz     " << MODz - zz[1] << endl
      ;

  } // MOD

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // done

  cout << endl << histoFile->GetName() << endl;

  cout << endl;

  return 0;
}
