
// Daniel Pitzl, DESY, Jun 2019
// telescope 3-D event display using ROOT

// needs runs.dat

// make ed3d
// ed3d 36663
// ed3d 36781 vertices

#include "eudaq/FileReader.hh"
#include "eudaq/PluginManager.hh"

#include <TApplication.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TH3.h>
#include <TPolyLine3D.h>
#include <TPolyMarker3D.h>

#include <sstream> // stringstream
#include <fstream> // filestream
#include <set>
#include <cmath> // fabs
#include <unistd.h> // usleep
#include <sys/ioctl.h> // kbhit

using namespace std;
using namespace eudaq;

struct pixel {
  int col;
  int row;
  int adc;
  double q;
  int ord;
  bool big;
};

struct cluster {
  vector <pixel> vpix; // Armin Burgmeier: list
  int size;
  int ncol, nrow;
  double col, row;
  double charge;
  bool big;
};

struct triplet {
  double xm;
  double ym;
  double zm;
  double sx;
  double sy;
  double rxy;
  int ghost;
  vector <double> vx;
  vector <double> vy;
};

struct vertex {
  double x;
  double y;
  double z;
  unsigned i;
  unsigned j;
};

//------------------------------------------------------------------------------
bool kbhit()
{
  usleep(1000); // [us]
  int byteswaiting;
  ioctl( 0, FIONREAD, &byteswaiting );
  return byteswaiting > 0;
}

//------------------------------------------------------------------------------
vector <cluster> getClus( vector <pixel> pb, int fCluCut = 1 ) // 1 = no gap
{
  // returns clusters with local coordinates
  // decodePixels should have been called before to fill pixel buffer pb 
  // simple clusterization
  // cluster search radius fCluCut ( allows fCluCut-1 empty pixels)

  vector<cluster> vc;
  if( pb.size() == 0 ) return vc;

  vector <bool> gone( pb.size() ); // initialized to zero

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
        if( !gone[i] ) { // unused pixel
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
    c.big = 0;
    int minx = 999;
    int maxx = 0;
    int miny = 999;
    int maxy = 0;

    for( vector<pixel>::iterator p = c.vpix.begin();  p != c.vpix.end();  ++p ) {
      double Qpix = p->q; // calibrated [Vcal]

      c.charge += Qpix;
      sumQ += Qpix;
      c.col += (*p).col*Qpix;
      c.row += (*p).row*Qpix;
      if( p->big ) c.big = 1;
      if( p->col > maxx ) maxx = p->col;
      if( p->col < minx ) minx = p->col;
      if( p->row > maxy ) maxy = p->row;
      if( p->row < miny ) miny = p->row;
    }

    //cout << "(cluster with " << c.vpix.size() << " pixels)" << endl;

    if( ! ( c.charge == 0 ) ) {
      c.col /= sumQ;
      c.row /= sumQ;
    }
    else {
      c.col = (*c.vpix.begin()).col;
      c.row = (*c.vpix.begin()).row;
      cout << "GetClus: cluster with zero charge" << endl;
    }

    c.ncol = maxx-minx+1;
    c.nrow = maxy-miny+1;

    vc.push_back(c); // add cluster to vector

    // look for a new seed = used pixel:

    while( ( ++seed < pb.size() ) && gone[seed] );

  } // while over seeds

  // nothing left, return clusters

  return vc;

} // getClus

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
  double pbeam = 5.6;

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

      stringstream tokenizer( line );
      string tag;
      tokenizer >> tag; // leading white space is suppressed
      if( tag.substr(0,1) == hash ) // comments start with #
	continue;

      if( tag == plane ) {
	tokenizer >> ipl;
	continue;
      }

      if( ipl < 0 || ipl >= 9 ) {
	cout << "wrong plane number " << ipl << endl;
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
      midx[ipl] = 0.5 * sizex[ipl]; // mid plane
      midy[ipl] = 0.5 * sizey[ipl]; // mid plane
    }

  } // geo scope

  geoFile.close();

  double X0m = 0.9E-3; // Mimosa and air
  double scatm = 0.0136 * sqrt(X0m) / pbeam * ( 1 + 0.038*log(X0m) ); // [rad] Mimosa scattering

  double kCut = 12 * scatm; // [rad]

  cout << endl << "kinkCut " << kCut*1E3 << " mrad" << endl;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // alignment:

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

      if( tag == iteration ) 
	tokenizer >> aligniteration;

      if( tag == plane )
	tokenizer >> ipl;

      if( ipl < 0 || ipl >= 9 ) {
	cout << "wrong plane number " << ipl << endl;
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

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // hot pixels:

  ostringstream hotFileName; // output string stream

  hotFileName << "hot_" << run << ".dat";

  ifstream ihotFile( hotFileName.str() );

  set <int> hotset[9];

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

      stringstream tokenizer( line );
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

  for( int ipl = 1; ipl <= 6; ++ipl )
    cout << ipl << ": hot " << hotset[ipl].size() << endl;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // ROOT:

  gStyle->SetTextFont( 62 ); // 62 = Helvetica bold
  gStyle->SetTextAlign( 11 );

  gStyle->SetTickLength( -0.02, "x" ); // tick marks outside
  gStyle->SetTickLength( -0.01, "y" );
  gStyle->SetTickLength( -0.01, "z" );

  gStyle->SetLabelOffset( 0.038, "x" );
  gStyle->SetLabelOffset( 0.022, "y" );
  gStyle->SetLabelOffset( 0.022, "z" );

  gStyle->SetTitleOffset( 1.2, "x" );
  gStyle->SetTitleOffset( 1.6, "y" );
  gStyle->SetTitleOffset( 1.1, "z" );

  gStyle->SetLabelFont( 62, "X" );
  gStyle->SetLabelFont( 62, "Y" );
  gStyle->SetLabelFont( 62, "z" );

  gStyle->SetTitleFont( 62, "X" );
  gStyle->SetTitleFont( 62, "Y" );
  gStyle->SetTitleFont( 62, "z" );

  gStyle->SetTitleBorderSize( 0 ); // no frame around global title
  gStyle->SetTitleAlign( 13 ); // 13 = left top align
  gStyle->SetTitleX( 0.12 ); // global title
  gStyle->SetTitleY( 0.98 ); // global title

  gStyle->SetLineWidth( 1 ); // frames
  gStyle->SetHistLineColor( 4 ); // 4=blau
  gStyle->SetHistLineWidth( 3 );
  gStyle->SetHistFillColor( 5 ); // 5=gelb
  //  gStyle->SetHistFillStyle(4050); // 4050 = half transparent
  gStyle->SetHistFillStyle( 1001 ); // 1001 = solid

  gStyle->SetFrameLineWidth( 2 );

  // statistics box:

  gStyle->SetOptStat( 0 );
  gStyle->SetStatFormat( "8.6g" ); // more digits, default is 6.4g
  gStyle->SetStatFont( 42 ); // 42 = Helvetica normal
  //  gStyle->SetStatFont(62); // 62 = Helvetica bold
  gStyle->SetStatBorderSize( 1 ); // no 'shadow'

  gStyle->SetStatX( 0.80 );
  gStyle->SetStatY( 0.95 );

  gStyle->SetPalette( 1 ); // rainbow colors

  gStyle->SetHistMinimumZero(  ); // no zero suppression

  gStyle->SetOptDate( 0 );

  TApplication app( "app", 0, 0 );

  TCanvas * c1 = new TCanvas( "c1", "RD53A event display", 500, 0, 900, 900 );

  c1->SetBottomMargin( 0.08 );
  c1->SetLeftMargin( 0.08 );
  c1->SetRightMargin( 0.03 );

  //gPad->Update();

  double dz = zz[6] - zz[1];

  TH3I xyzview( "",
		"display;x [mm];z [mm];y [mm]",
		110, -11, 11, 120, -0.1*dz, 1.1*dz, 120, -6, 6 );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // event loop:

  cout << endl;

  int iev = 0;

  int nevd = 0;
  uint64_t evTLU0 = 0;
  const double fTLU = 384E6; // 384 MHz TLU clock

  string MIM{"MIMOSA26"};

  bool more = 1; // event displays
  string Q{"q"};

  do {
    // Get next event:
    DetectorEvent evt = reader->GetDetectorEvent();

    if( evt.IsBORE() ) {
      eudaq::PluginManager::Initialize(evt);
      continue;
    }

    ++iev;

    bool ldbg = 0;
    if( lev < 99 )
      ldbg = 1;

    uint64_t evTLU = evt.GetTimestamp(); // 384 MHz = 2.6 ns
    if( iev < 2  )
      evTLU0 = evTLU;
    double evsec = (evTLU - evTLU0) / fTLU;

    if( iev < 10 || ldbg )
      cout << "ed3d processing  " << iev << "  taken " << evsec << endl;
    else if( iev < 100 && iev%10 == 0 )
      cout << "ed3d processing  " << iev << "  taken " << evsec << endl;
    else if( iev < 1000 && iev%100 == 0 )
      cout << "ed3d processing  " << iev << "  taken " << evsec << endl;
    else if( iev%1000 == 0 )
      cout << "ed3d processing  " << iev << "  taken " << evsec << endl;

    StandardEvent sevt = eudaq::PluginManager::ConvertToStandard(evt);

    vector <cluster> cl[9];
    int ncl = 0;

    int mpl = 1; // plane numbering start at 1

    for( size_t ipl = 0; ipl < sevt.NumPlanes(); ++ipl ) {

      const eudaq::StandardPlane &plane = sevt.GetPlane(ipl);

      std::vector<double> pxl = plane.GetPixels<double>();

      if( ldbg )
	cout << "  " << ipl
	     << ": plane " << plane.ID()
	     << " " << plane.Type() // NI
	     << " " << plane.Sensor() // MIMOSA26
	     << " frames " << plane.NumFrames() // 2
	     << " pivot " << plane.PivotPixel() // 6830
	     << " total " << plane.TotalPixels() // 663552
	     << " hits " << plane.HitPixels() // 5
	     << flush;

      if( plane.Sensor() != MIM ) {
	if( ldbg ) cout << endl;
	continue;
      }

      vector <pixel> pb; // for clustering

      for( size_t ipix = 0; ipix < pxl.size(); ++ipix ) {

	if( ldbg ) 
	  cout << plane.GetX(ipix)
	       << " " << plane.GetY(ipix)
	       << " " << plane.GetPixel(ipix) << " "
	       << flush;

	int ix = plane.GetX(ipix); // global column 0..415
	int iy = plane.GetY(ipix); // global row 0..159
	int adc = plane.GetPixel(ipix); // ADC 0..255

	// skip hot pixels:

	int ipx = ix*ny[mpl] + iy;
	if( hotset[mpl].count(ipx) ) continue;

	double q = adc;

	// fill pixel block for clustering:

	pixel px;
	px.col = ix; // col
	px.row = iy; // row
	px.adc = adc;
	px.q = q;
	px.ord = pb.size(); // readout order
	px.big = 0;
	pb.push_back(px);

      } // pix

      if( ldbg ) cout << endl;

      // clustering:

      cl[mpl] = getClus( pb );

      if( ldbg ) cout << "  z " << int(zz[mpl]+0.5)
		      << ", clusters " << cl[mpl].size()
		      << endl;

      ncl += cl[mpl].size();

      ++mpl;
      if( mpl > 6 ) break; // skip others

    } // planes

    // event selection:

    if( ncl < 6 ) continue;

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // make triplets 2+0-1:

    vector <triplet> triplets;

    double z1 = zz[1] + alignz[1];
    double z2 = zz[2] + alignz[2];
    double z3 = zz[3] + alignz[3];
    double dz13 = z3 - z1;

    for( vector<cluster>::iterator cA = cl[1].begin(); cA != cl[1].end(); ++cA ) {

      double xA = cA->col*ptchx[1] - alignx[1];
      double yA = cA->row*ptchy[1] - aligny[1];
      double xmid = xA - midx[1];
      double ymid = yA - midy[1];
      xA = xmid - ymid*rotx[1];
      yA = ymid + xmid*roty[1];

      for( vector<cluster>::iterator cC = cl[3].begin(); cC != cl[3].end(); ++cC ) {

	double xC = cC->col*ptchx[3] - alignx[3];
	double yC = cC->row*ptchy[3] - aligny[3];
	double xmid = xC - midx[3];
	double ymid = yC - midy[3];
	xC = xmid - ymid*rotx[3];
	yC = ymid + xmid*roty[3];
 
	double slpx = ( xC - xA ) / dz13; // slope x
	double slpy = ( yC - yA ) / dz13; // slope y

	double avx = 0.5 * ( xA + xC ); // mid
	double avy = 0.5 * ( yA + yC );
	double avz = 0.5 * ( z1 + z3 ); // mid z

	// middle plane B = 2:

	for( vector<cluster>::iterator cB = cl[2].begin(); cB != cl[2].end(); ++cB ) {

	  double xB = cB->col*ptchx[2] - alignx[2];
	  double yB = cB->row*ptchy[2] - aligny[2];
	  double xmid = xB - midx[2];
	  double ymid = yB - midy[2];
	  xB = xmid - ymid*rotx[2];
	  yB = ymid + xmid*roty[2];

	  // kinks:

	  double kx = (xC-xB)/(z3-z2) - (xB-xA)/(z2-z1);
	  double ky = (yC-yB)/(z3-z2) - (yB-yA)/(z2-z1);

	  if( fabs( kx ) > kCut ) continue;
	  if( fabs( ky ) > kCut ) continue;

          // interpolate track to B:

          double dz = z2 - avz;
          double xk = avx + slpx * dz; // triplet at k
          double yk = avy + slpy * dz;

          double dx3 = xB - xk;
          double dy3 = yB - yk;

	  triplet tri;
	  tri.xm = avx;
	  tri.ym = avy;
	  tri.zm = avz;
	  tri.sx = slpx;
	  tri.sy = slpy;
	  tri.rxy = sqrt(dx3*dx3+dy3*dy3);
	  tri.ghost = 0;

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

	  triplets.push_back(tri);

	} // cl B

      } // cl C

    } // cl A

    if( ldbg ) cout << "triplets " << triplets.size() << endl;

    // event selection:

    if( triplets.size() < 1 ) continue;

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // make driplets 3+5-4:

    vector <triplet> driplets;

    double z4 = zz[4] + alignz[4];
    double z5 = zz[5] + alignz[5];
    double z6 = zz[6] + alignz[6];
    double dz46 = z6 - z4;

    for( vector<cluster>::iterator cA = cl[4].begin(); cA != cl[4].end(); ++cA ) {

      double xA = cA->col*ptchx[4] - alignx[4];
      double yA = cA->row*ptchy[4] - aligny[4];
      double xmid = xA - midx[4];
      double ymid = yA - midy[4];
      xA = xmid - ymid*rotx[4];
      yA = ymid + xmid*roty[4];

      for( vector<cluster>::iterator cC = cl[6].begin(); cC != cl[6].end(); ++cC ) {

	double xC = cC->col*ptchx[6] - alignx[6];
	double yC = cC->row*ptchy[6] - aligny[6];
	double xmid = xC - midx[6];
	double ymid = yC - midy[6];
	xC = xmid - ymid*rotx[6];
	yC = ymid + xmid*roty[6];
 
	double slpx = ( xC - xA ) / dz46; // slope x
	double slpy = ( yC - yA ) / dz46; // slope y

	double avx = 0.5 * ( xA + xC ); // mid
	double avy = 0.5 * ( yA + yC );
	double avz = 0.5 * ( z4 + z6 ); // mid z

	// middle plane B = 5:

	for( vector<cluster>::iterator cB = cl[5].begin(); cB != cl[5].end(); ++cB ) {

	  double xB = cB->col*ptchx[5] - alignx[5];
	  double yB = cB->row*ptchy[5] - aligny[5];
	  double xmid = xB - midx[5];
	  double ymid = yB - midy[5];
	  xB = xmid - ymid*rotx[5];
	  yB = ymid + xmid*roty[5];

	  // kinks:

	  double kx = (xC-xB)/(z6-z5) - (xB-xA)/(z5-z4);
	  double ky = (yC-yB)/(z6-z5) - (yB-yA)/(z5-z4);

	  if( fabs( kx ) > kCut ) continue;
	  if( fabs( ky ) > kCut ) continue;

          // interpolate track to B:

          double dz = z5 - avz;
          double xk = avx + slpx * dz; // triplet at k
          double yk = avy + slpy * dz;

          double dx3 = xB - xk;
          double dy3 = yB - yk;

	  triplet dri;
	  dri.xm = avx;
	  dri.ym = avy;
	  dri.zm = avz;
	  dri.sx = slpx;
	  dri.sy = slpy;
	  dri.rxy = sqrt(dx3*dx3+dy3*dy3);
	  dri.ghost = 0;

	  vector <double> ux(3);
	  ux[0] = xA;
	  ux[1] = xB;
	  ux[2] = xC;
	  dri.vx = ux;

	  vector <double> uy(3);
	  uy[0] = yA;
	  uy[1] = yB;
	  uy[2] = yC;
	  dri.vy = uy;

	  driplets.push_back(dri);

	} // cl B

      } // cl C

    } // cl A

    if( ldbg ) cout << "driplets " << driplets.size() << endl;

    // event selection:

    if( driplets.size() < 1 ) continue;

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

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // driplet vertices:

    vector <vertex> vvt;
    vector <vertex> vv0;

    for( unsigned int jB = 0; jB < driplets.size(); ++jB ) { // j = B = downstream

      if( driplets[jB].ghost ) continue;

      double xmB = driplets[jB].xm;
      double ymB = driplets[jB].ym;
      double zmB = driplets[jB].zm;
      double sxB = driplets[jB].sx;
      double syB = driplets[jB].sy;

      // vertices:

      for( unsigned int jj = jB + 1; jj < driplets.size(); ++jj ) {

	if( driplets[jj].ghost ) continue;

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

	double dd = ( (xmj-xmB)*Nx + (ymj-ymB)*Ny ) / Nn; // signed

	if( fabs(dd) > 0.040 ) continue;

	// z intersection:

	double oax = sxj - sxB;
	double oay = syj - syB;
	double zix = ( xmB - xmj ) / oax;
	double ziy = ( ymB - ymj ) / oay;

	// weight by opening angle:

	double zv = zmB + ( oax*oax*zix + oay*oay*ziy ) / (oax*oax+oay*oay);

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

	  if( fabs(dxv) < 0.060 && fabs(dyv) < 0.060 )
	    ++nmt;

	} // iA

	vertex vx;
	vx.x = xv;
	vx.y = yv;
	vx.z = zv;
	vx.i = jB;
	vx.j = jj;
	if( nmt )
	  vvt.push_back( vx );
	else
	  vv0.push_back( vx );

      } // jj

    } // jB

    // event selection:

    //if( vvt.size() + vv0.size() == 0 ) continue;

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // display

    ++nevd;

    cout
      << "event " << iev
      << " is display " << nevd
      << " with " << ncl << " clusters"
      << ", " << triplets.size() << " triplets"
      << ", " << driplets.size() << " driplets"
      << ", " << vvt.size()+vv0.size() << " vertices"
      << endl;

    xyzview.SetTitle( Form( "event %i", iev ) );

    xyzview.Draw();

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // vertices:

    double xyzt[3*vvt.size()];

    for( unsigned int iv = 0; iv < vvt.size(); ++iv ) {
      xyzt[3*iv+0] = vvt[iv].x;
      xyzt[3*iv+1] = vvt[iv].z;
      xyzt[3*iv+2] = vvt[iv].y;
      cout << "  vx"
	   << " " << vvt[iv].i
	   << "-" << vvt[iv].j
	   << " " << vvt[iv].x
	   << " " << vvt[iv].y
	   << " " << vvt[iv].z
	   << " with track"
	   << endl;
    }

    TPolyMarker3D * pmvt = new TPolyMarker3D( vvt.size(), xyzt );
    pmvt->SetMarkerColor(1);
    pmvt->SetMarkerStyle(25);
    pmvt->SetMarkerSize(1.5);
    pmvt->Draw(); // without axis option: overlay

    double xyz0[3*vv0.size()];

    for( unsigned int iv = 0; iv < vv0.size(); ++iv ) {
      xyz0[3*iv+0] = vv0[iv].x;
      xyz0[3*iv+1] = vv0[iv].z;
      xyz0[3*iv+2] = vv0[iv].y;
      cout << "  vx"
	   << " " << vv0[iv].i
	   << "-" << vv0[iv].j
	   << " " << vv0[iv].x
	   << " " << vv0[iv].y
	   << " " << vv0[iv].z
	   << " without"
	   << endl;
    }

    TPolyMarker3D * pmv0 = new TPolyMarker3D( vv0.size(), xyz0 );
    pmv0->SetMarkerColor( 8);
    pmv0->SetMarkerStyle(25);
    pmv0->SetMarkerSize(1.5);
    pmv0->Draw(); // without axis option: overlay

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // triplets:

    vector <TPolyLine3D> lines;

    double lx[2];
    double ly[2];
    double lz[2];

    double zb = z1 - 0.5 * ( z2 - z1 ); // begin track
    double zm = 0.5 * ( z3 + z4 ); // mid

    for( unsigned int iA = 0; iA < triplets.size(); ++iA ) { // iA = upstream

      if( triplets[iA].ghost ) continue;

      double avx = triplets[iA].xm;
      double avy = triplets[iA].ym;
      double avz = triplets[iA].zm;
      double slx = triplets[iA].sx;
      double sly = triplets[iA].sy;

      double d1 = zb - avz;
      double x1 = avx + slx * d1;
      double y1 = avy + sly * d1;
      lx[0] = x1;
      ly[0] = zb;
      lz[0] = y1;

      double d3 = zm - avz;
      double x3 = avx + slx * d3;
      double y3 = avy + sly * d3;
      lx[1] = x3;
      ly[1] = zm;
      lz[1] = y3;

      TPolyLine3D ll( 2, lx, ly, lz );
      ll.SetLineColor(6);
      ll.SetLineWidth(2);
      lines.push_back(ll);

    } // loop triplets iA

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // driplets:

    double ze = z6 + 0.5 * ( z6 - z5 ); // end

    for( unsigned int jB = 0; jB < driplets.size(); ++jB ) { // jB = downstream

      if( driplets[jB].ghost ) continue;

      double avx = driplets[jB].xm;
      double avy = driplets[jB].ym;
      double avz = driplets[jB].zm;
      double slx = driplets[jB].sx;
      double sly = driplets[jB].sy;

      double d4 = zm - avz;
      double x4 = avx + slx * d4;
      double y4 = avy + sly * d4;
      lx[0] = x4;
      ly[0] = zm;
      lz[0] = y4;

      double d6 = ze - avz;
      double x6 = avx + slx * d6;
      double y6 = avy + sly * d6;
      lx[1] = x6;
      ly[1] = ze;
      lz[1] = y6;

      TPolyLine3D ll( 2, lx, ly, lz );
      ll.SetLineColor(7);
      ll.SetLineWidth(2);
      lines.push_back(ll);

    } // loop driplets

    for( size_t ii = 0; ii < lines.size(); ++ii )
      lines[ii].Draw( "same" );

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Mimosa hits:

    TPolyMarker3D * pm3d[6]; // per plane

    for( int ipl = 1; ipl <= 6; ++ipl ) {

      double hits[3*cl[ipl].size()]; // xyz clusters, all planes

      int icl = 0;

      for( vector<cluster>::iterator cA = cl[ipl].begin(); cA != cl[ipl].end(); ++cA ) {

	double xA = cA->col*ptchx[ipl] - alignx[ipl];
	double yA = cA->row*ptchy[ipl] - aligny[ipl];

	double zA = zz[ipl] + alignz[ipl];

	double xmid = xA - midx[ipl];
	double ymid = yA - midy[ipl];

	xA = xmid - ymid*rotx[ipl];
	yA = ymid + xmid*roty[ipl];

	hits[3*icl+0] = xA;
	hits[3*icl+1] = zA;
	hits[3*icl+2] = yA;

	++icl;

      } // cl

      pm3d[ipl-1] = new TPolyMarker3D( icl, hits );
      if(      ipl == 1 )
	pm3d[ipl-1]->SetMarkerColor(99);
      else if( ipl == 2 )
	pm3d[ipl-1]->SetMarkerColor(96);
      else if( ipl == 3 )
	pm3d[ipl-1]->SetMarkerColor(93);
      else if( ipl == 4 )
	pm3d[ipl-1]->SetMarkerColor( 4);
      else if( ipl == 5 )
	pm3d[ipl-1]->SetMarkerColor(51);
      else
	pm3d[ipl-1]->SetMarkerColor(64);
      pm3d[ipl-1]->SetMarkerStyle(20);
      pm3d[ipl-1]->SetMarkerSize(0.8);
      pm3d[ipl-1]->Draw(); // without axis option: overlay

    } // planes

    c1->Update();

    cout << "enter any key, q to stop" << endl;

    while( !kbhit() )
      gSystem->ProcessEvents(); // ROOT

    string any;
    cin >> any;
    if( any == Q )
      more = 0;

  } while( more && iev < lev && reader->NextEvent() );

  delete reader;

  cout << "done after " << iev << " events" << endl;

  cout << endl;

  return 0;
}
