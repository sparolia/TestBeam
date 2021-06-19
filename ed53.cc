
// Daniel Pitzl, DESY, Oct 2018
// event display for RD53A from eudaq

// make ed53
// ed53 34239
// ed53 36619
// needs runs.dat

#include "eudaq/FileReader.hh"
#include "eudaq/PluginManager.hh"
#include "../main/lib/plugins/BDAQ53ConverterPlugin.cc"

#include <TApplication.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TH2.h>

#include <map>
#include <sys/ioctl.h>

using namespace std;
using namespace eudaq;

struct pixel {
  int col;
  int row;
  int tot;
  int frm;
  bool pivot;
};

struct cluster {
  vector <pixel> vpix; // Armin Burgmeier: list
  int size;
  int ncol, nrow, nfrm;
  double col, row;
  int signal;
  double mindxy;
};

//------------------------------------------------------------------------------
bool kbhit()
{
  int byteswaiting;
  ioctl( 0, FIONREAD, &byteswaiting );
  return byteswaiting > 0;
}

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
}

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

  int lev = 999222111; // last event

  for( int i = 1; i < argc; ++i ) {

    if( !strcmp( argv[i], "-l" ) )
      lev = atoi( argv[++i] ); // last event

  } // argc

  gStyle->SetTextFont(62); // 62 = Helvetica bold
  gStyle->SetTextAlign(11);

  gStyle->SetTitleBorderSize(0); // no frame around global title
  gStyle->SetTitleAlign(13); // 13 = left top align
  gStyle->SetTitleX( 0.15 ); // global title
  gStyle->SetTitleY( 0.98 ); // global title

  gStyle->SetTitleFont( 62, "XYZ" );

  gStyle->SetTitleOffset( 1.3, "x" );
  gStyle->SetTitleOffset( 1.5, "y" );
  gStyle->SetTitleOffset( 1.5, "z" );

  gStyle->SetLabelFont( 62, "XYZ" );

  gStyle->SetLabelOffset( 0.022, "xyz" );

  gStyle->SetTickLength( -0.02, "xyz" ); // tick marks outside

  gStyle->SetLineWidth(1);// frames
  gStyle->SetHistLineColor(4); // 4=blau
  gStyle->SetHistLineWidth(3);
  gStyle->SetHistFillColor(5); // 5=gelb
  //  gStyle->SetHistFillStyle(4050); // 4050 = half transparent
  gStyle->SetHistFillStyle(1001); // 1001 = solid

  gStyle->SetFrameLineWidth(2);

  // statistics box:

  gStyle->SetOptStat(10);
  gStyle->SetStatFont(42); // 42 = Helvetica normal
  gStyle->SetStatBorderSize(1); // no 'shadow'
  gStyle->SetStatX(0.82);
  gStyle->SetStatY(0.92);

  gStyle->SetPalette(55); // sunset
  //gStyle->SetPalette(56); // white to red
  //gStyle->SetPalette(73); // blue
  //gStyle->SetPalette(90); // green to magenta
  //gStyle->SetPalette(109); // sea blue to magenta
  gStyle->SetNumberContours(16); // 0..16

  TApplication app( "app", 0, 0 );
  TCanvas c1( "c1", "RD53A event display", 900, 800 ); // square
  c1.SetTopMargin( 0.12 );
  c1.SetRightMargin( 0.18 );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // runs.dat:

  cout << endl;

  int chip0 = 501;
  int fifty = 1; // default is 50x50
  int rot90 = 0; // default is straight

  ifstream runsFile( "runs.dat" );

  if( runsFile.bad() || ! runsFile.is_open() ) {
    cout << "Error opening runs.dat" << endl;
    return 1;
  }
  // can there be instructions between if and else ? no

  else {

    cout << "read runs from runs.dat" << endl;

    string hash( "#" );
    string RUN( "run" );
    string CHIP( "chip" );
    string FIFTY( "fifty" );
    string ROT90( "rot90" );
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

      if( tag == CHIP ) {
	tokenizer >> chip0;
	continue;
      }

      if( tag == FIFTY ) {
	tokenizer >> fifty;
	continue;
      }

      if( tag == ROT90 ) {
	tokenizer >> rot90;
	continue;
      }

      // anything else on the line and in the file gets ignored

    } // while getline

    if( found )
      cout 
	<< "settings for run " << run << ":" << endl
	<< "  DUT chip " << chip0 << endl
	<< "  fifty " << fifty << endl
	<< "  rot90 " << rot90 << endl
	<< endl;
    else {
      cout << "run " << run << " not found in runs.dat" << endl;
      return 1;
    }

  } // runsFile

  runsFile.close();

  int nbc = 400;
  int nbr = 192;
  if( !fifty ) {
    nbc = 200;
    nbr = 384;
  }

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

  int iev = 0;

  bool ldbg = 0;
  bool more = 1;

  do {

    DetectorEvent evt = reader->GetDetectorEvent();

    if( evt.IsBORE() )
      eudaq::PluginManager::Initialize(evt);

    StandardEvent sevt = eudaq::PluginManager::ConvertToStandard(evt);

    if( ldbg ) cout << "planes " << sevt.NumPlanes() << endl;

    vector < cluster > vcl;

    map <int, int> cols;
    map <int, int> rows;

    for( size_t iplane = 0; iplane < sevt.NumPlanes(); ++iplane ) {

      const eudaq::StandardPlane &plane = sevt.GetPlane(iplane);

      if( ldbg )
	std::cout
	  << "  " << iplane
	  << ": plane " << plane.ID() // 1
	  << " " << plane.Type() // NI or BDAQ53
	  << " " << plane.Sensor() // MIMOSA26 or RD53A
	  << " frames " << plane.NumFrames() // 2 for NI or 32 for BDAQ53
	  << " pivot " << plane.PivotPixel() // 4166 or 0
	  << " total " << plane.TotalPixels()
	  << " hits " << plane.HitPixels()
	  ;

      int ipl = plane.ID(); // 0 = RD53A

      if( ipl > 0 )
	continue;

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

	  pixel px;
	  px.col = ix; // col
	  px.row = iy; // row
	  px.tot = tot + 1; // shift from zero
	  px.frm = frm;
	  px.pivot = plane.GetPivot(ipix,frm);

	  // map from ROC to sensor:

	  if( fifty ) {
	    px.col = ix; // 50x50
	    px.row = iy;
	  }
	  else {
	    px.col = ix/2; // 100 um
	    if( ix%2 )
	      px.row = 2*iy + 0; // different from R4S
	    else
	      px.row = 2*iy + 1;
	  }

	  pb.push_back(px);

	  ++cols[px.col];
	  ++rows[px.row];

	} // pix

      if( ldbg )
	cout << endl;

      // clustering:

      vcl = getClusq( pb );

      if( ldbg )
	cout << "    clusters " << vcl.size() << endl;

    } // planes

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // DUT:

    TH2D hpxmap( "pxmap",
		 Form( "pixel map %i;sensor column;sensor row;signal [ToT]", iev ),
		 nbc, -0.5, nbc-0.5, nbr, -0.5, nbr-0.5 );
    hpxmap.SetMinimum(0);
    hpxmap.SetMaximum(16);

    int maxncol = 0;
    int maxnrow = 0;
    for( vector<cluster>::iterator c = vcl.begin(); c != vcl.end(); ++c ) {

      if( c->ncol > maxncol )
	maxncol = c->ncol;

      if( c->nrow > maxnrow )
	maxnrow = c->nrow;

      for( vector<pixel>::iterator px = c->vpix.begin(); px != c->vpix.end(); ++px ) {

	hpxmap.Fill( px->col, px->row, px->tot );

      } // px

    } // clus

    //if( maxncol > 33 ) {
    if( maxncol > 22 ) {
      //if( rows.size() > 22 ) {

      hpxmap.Draw( "colz" );
      c1.Update();

      cout << "event " << iev
	   << ". enter any key, q to stop" << endl;

      while( !kbhit() ) // ioctl
	gSystem->ProcessEvents(); // ROOT

      string q {"q"};
      string any;
      cin >> any;
      if( any == q )
	more = 0;

    } // show

    ++iev;

  } while( more && reader->NextEvent() && iev < lev );

  delete reader;

  cout << "done after " << iev << " events" << endl;

  cout << endl;

  return 0;
}
