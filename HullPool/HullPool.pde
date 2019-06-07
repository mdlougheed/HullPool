// A JBOX2D buoyant body wavepool simulation
// by Mark D. Lougheed
// Project STRUIX - 2015


// JBOX2D Library Components
import shiffman.box2d.*;
import org.jbox2d.common.*;
import org.jbox2d.dynamics.joints.*;
import org.jbox2d.collision.shapes.*;
import org.jbox2d.collision.shapes.Shape;
import org.jbox2d.dynamics.*;
import org.jbox2d.dynamics.contacts.*;

// G4P GUI Library Components
import g4p_controls.*;
GCustomSlider sdr;



// A reference to our box2d world
Box2DProcessing box2d;
float t; // Time counter


// Objects in the world
WaveBody surface;
WaveTank tank;

ArrayList<BuoyantBody> theHull;
HullBody hull;


float dt;
Vec2 dv;
Vec2 prevV;


PrintWriter output;

PIDController PID;

//float maxThrust=11120.0;  // Maximum thrust in N
float maxThrust=25000.0;  // Maximum thrust in N

float targetMass= 14028.0;  //  Hull mass - kg
float targetVolume=72.6195339; // Total Hull Volume - m^3

float sigWaveHt=0.0;


void setup() {

    MassData mData;
  float polyArea;
  
  
  
  Settings.maxPolygonVertices=100;
  //Settings.maxPolygonVertices=10;  // User control of maximum PolygonShape vertices???





  size(1000,600);
  
  smooth();

  // Initialize G4P GUI
  makeSliderConfigControls();


  // Initialize box2d physics and create the world
  box2d = new Box2DProcessing(this);
  box2d.createWorld();
  

  
  
  // Set world to earth gravity
  box2d.setGravity(0, -9.80665);
  
  // Turn on collision listening!
  box2d.listenForCollisions();

  // Time starts at zero
  t=0;
  dt = 1f / 100f;

dv=new Vec2(0,0);
prevV = new Vec2(0,0);

// Add some stuff to the world

  // Create the WaveBody
  surface = new WaveBody(-40,0);
//  surface.SetWaveParams(8.2, 2.5);  // Wave period (s), H13 (m)
//  surface.SetWaveParams(5, 2);  // Wave period (s), H13 (m)  
  surface.SetWaveParams(6, 1.0);  // Wave period (s), H13 (m)
  surface.fetchLength=1000;

  
  
 // Create the tank
 // WaveTank(float posX, float posY, float _depth, float _length, float _borderThickness) {
  tank = new WaveTank(-41,-11, 12,1000, 1);
  

  // Create a vessel HullBody to contain all the BuoyantBodies (slices) that make up the hull  
  theHull = new ArrayList<BuoyantBody>();


  // Shape outlines in local coordinates
  // defined counter-clockwise (interior is to the left of each edge)
  ArrayList<Vec2> outline = new ArrayList<Vec2>();
  
// Slice #1
  outline.clear();

     outline.add(new Vec2(-0.0054, 0.0619)); // 0
     outline.add(new Vec2(0.0149, 0.0619)); // 1
     outline.add(new Vec2(1.6028, 0.0619)); // 2
     outline.add(new Vec2(2.6122, 0.0620)); // 3
     outline.add(new Vec2(3.2472, 0.0619)); // 4
     outline.add(new Vec2(4.3047, 0.0617)); // 5
     outline.add(new Vec2(4.9903, 0.0620)); // 6
     outline.add(new Vec2(6.4208, 0.0625)); // 7
     outline.add(new Vec2(6.7332937231109558, 0.061888119248493703)); // 8
     outline.add(new Vec2(7.0814381186832449, 0.060900899890078172)); // 9
     outline.add(new Vec2(7.4512022418962811, 0.061509691407730743)); // 10
     outline.add(new Vec2(7.636824622149633, 0.06319824186622347)); // 11
     outline.add(new Vec2(7.8228043764736386, 0.066218961310349558)); // 12
     outline.add(new Vec2(8.005525811852408, 0.070772162461973931)); // 13
     outline.add(new Vec2(8.0091696339144427, 0.070880891333214374)); // 14
     outline.add(new Vec2(8.1976403691122091, 0.077557454308143425)); // 15
     outline.add(new Vec2(8.3889672342158601, 0.08663214276312059)); // 16
     outline.add(new Vec2(8.5817924024689631, 0.09840255292068148)); // 17
     outline.add(new Vec2(8.7764686117132733, 0.11330604841294026)); // 18
     outline.add(new Vec2(8.9726246308163553, 0.13195354008279817)); // 19
     outline.add(new Vec2(9.100595202927666, 0.14639656312537955)); // 20
     outline.add(new Vec2(9.1681463960734177, 0.15482709552979562)); // 21
     outline.add(new Vec2(9.3618518570655702, 0.1824264677501862)); // 22
     outline.add(new Vec2(9.5528448506442167, 0.21522655919471984)); // 23
     outline.add(new Vec2(9.7404230434961381, 0.25367887560585933)); // 24
     outline.add(new Vec2(9.9240589500732206, 0.29822057398989854)); // 25
     outline.add(new Vec2(10.078839775711256, 0.34181817307712092)); // 26
     outline.add(new Vec2(10.103335577961053, 0.34927473692540295)); // 27
     outline.add(new Vec2(10.279927725898835, 0.40793349568028336)); // 28
     outline.add(new Vec2(10.455142487597897, 0.47541295440316406)); // 29
     outline.add(new Vec2(10.627050476348625, 0.55153952719215416)); // 30
     outline.add(new Vec2(10.793578547628133, 0.63490606946921635)); // 31
     outline.add(new Vec2(10.954757630735571, 0.72468585500179339)); // 32
     outline.add(new Vec2(11.021052844041918, 0.76422696566201365)); // 33
     outline.add(new Vec2(11.110503264271095, 0.82000715635163612)); // 34
     outline.add(new Vec2(11.261556532780773, 0.92054931712456378)); // 35
     outline.add(new Vec2(11.408108251362375, 1.0257263658679321)); // 36
     outline.add(new Vec2(11.55019035810642, 1.1348744646877442)); // 37
     outline.add(new Vec2(11.687777365618075, 1.247299818813131)); // 38
     outline.add(new Vec2(11.774939647139318, 1.3219428061924905)); // 39
     outline.add(new Vec2(12.268955141902158, 1.7310922339513186)); // 40
     outline.add(new Vec2(12.7468, 2.1589)); // 41
     outline.add(new Vec2(-0.1885, 2.1589)); // 42
  
  float beam=0.4f; // 2x slice width
  Vec2 boatCG = new Vec2(4.59, 0.807);

  theHull.add(new BuoyantBody(100, 375, outline, beam));

  // suss out polygon area by using a MassData object
  // by using a density of 1.0 - the computed "mass" will equal the polygon area.
  mData = new MassData();  
  ((PolygonShape)theHull.get(theHull.size()-1).body.getFixtureList().getShape()).computeMass(mData, 1.0); 
  polyArea = mData.mass;
  
  // Slice mass is its fractional contribution of total hull volume
  theHull.get(theHull.size()-1).setBuoyantBodyMass( ((polyArea*beam)/targetVolume) * targetMass ); // CG will be reset to center of area by default.
  theHull.get(theHull.size()-1).setBuoyantBodyCG(boatCG); // be sure to set CG after setting the mass

  // Propulsor point.
  // 0,0 = through CG
  theHull.get(theHull.size()-1).setBuoyantBodyPropulsorPoint(/* boatCG */ new Vec2(-0.5, -0.5));


 
// Slice #2
  outline.clear();
  
     outline.add(new Vec2(-0.010819004580252766, 0.12389505242952875)); // 0
     outline.add(new Vec2(0.029874317929242257, 0.12390253790635553)); // 1
     outline.add(new Vec2(1.4285805788374972, 0.12390432553878133)); // 2
     outline.add(new Vec2(2.3177406541576966, 0.1239535473631023)); // 3
     outline.add(new Vec2(2.9403260527661383, 0.12389472771420718)); // 4
     outline.add(new Vec2(3.9772781828177695, 0.12373628047274313)); // 5
     outline.add(new Vec2(4.6495307644347026, 0.12402368633661377)); // 6
     outline.add(new Vec2(6.0520714184401303, 0.12440866215943958)); // 7
     outline.add(new Vec2(6.3583267877526666, 0.12370365927691603)); // 8
     outline.add(new Vec2(6.6995748210106392, 0.12274015004663952)); // 9
     outline.add(new Vec2(7.0841019068588453, 0.12337495856508957)); // 10
     outline.add(new Vec2(7.2778089850327863, 0.12513721003288603)); // 11
     outline.add(new Vec2(7.4721783694337063, 0.12829438564371715)); // 12
     outline.add(new Vec2(7.663430591981383, 0.13306040395109503)); // 13
     outline.add(new Vec2(7.6672455377400368, 0.13317423930800107)); // 14
     outline.add(new Vec2(7.8664296588182818, 0.14023328649659558)); // 15
     outline.add(new Vec2(8.0712112565524379, 0.14994818298967502)); // 16
     outline.add(new Vec2(8.2788538970158072, 0.16262379195651805)); // 17
     outline.add(new Vec2(8.4901009938281007, 0.17879934705717757)); // 18
     outline.add(new Vec2(8.7044828960609379, 0.19918051905975542)); // 19
     outline.add(new Vec2(8.8444008800974991, 0.21497187632659076)); // 20
     outline.add(new Vec2(8.9181030059912452, 0.22417012791724944)); // 21
     outline.add(new Vec2(9.1289287524677718, 0.2542063643760839)); // 22
     outline.add(new Vec2(9.3354993388798437, 0.28967771809167608)); // 23
     outline.add(new Vec2(9.5367416282875546, 0.33092588530286826)); // 24
     outline.add(new Vec2(9.7319321792812019, 0.37826386477616553)); // 25
     outline.add(new Vec2(9.8948296066028618, 0.42414171837246617)); // 26
     outline.add(new Vec2(9.920568536718477, 0.43197650642261598)); // 27
     outline.add(new Vec2(10.106331827764315, 0.49368675753195823)); // 28
     outline.add(new Vec2(10.2921677695279, 0.56526551115535761)); // 29
     outline.add(new Vec2(10.474501962293862, 0.64599952793667526)); // 30
     outline.add(new Vec2(10.649283866789691, 0.73348944194239085)); // 31
     outline.add(new Vec2(10.816623935404277, 0.82669197162366759)); // 32
     outline.add(new Vec2(10.884637589797411, 0.86725494746945397)); // 33
     outline.add(new Vec2(10.976401839128524, 0.92447383302028618)); // 34
     outline.add(new Vec2(11.130148340845157, 1.0268038842601614)); // 35
     outline.add(new Vec2(11.278293664540316, 1.1331203533542851)); // 36
     outline.add(new Vec2(11.420950280185629, 1.2427057354309401)); // 37
     outline.add(new Vec2(11.451678262959586, 1.2672402516334391)); // 38
     outline.add(new Vec2(12.043936356166103, 1.7002329661365794)); // 39
     outline.add(new Vec2(12.616461594405912, 2.1589999999999998)); // 40
     outline.add(new Vec2(-0.1885323943594, 2.1589999999999989)); // 41
  
  theHull.add(new BuoyantBody(300, 300, outline, beam));

  // suss out polygon area by using a MassData object
  // by using a density of 1.0 - the computed "mass" will equal the polygon area.
  mData = new MassData();  
  ((PolygonShape)theHull.get(theHull.size()-1).body.getFixtureList().getShape()).computeMass(mData, 1.0); 
  polyArea = mData.mass;
  
  // Slice mass is its fractional contribution of total hull volume
  theHull.get(theHull.size()-1).setBuoyantBodyMass( ((polyArea*beam)/targetVolume) * targetMass ); // CG will be reset to center of area by default.
  theHull.get(theHull.size()-1).setBuoyantBodyCG(boatCG); // be sure to set CG after setting the mass


// Slice #3
  outline.clear();
  
     outline.add(new Vec2(-0.016228506870379258, 0.18584257864429304)); // 0
     outline.add(new Vec2(0.044811476893863444, 0.18585380685953329)); // 1
     outline.add(new Vec2(1.2543382880922569, 0.18585648830817209)); // 2
     outline.add(new Vec2(2.023235755951541, 0.18589864987913948)); // 3
     outline.add(new Vec2(2.6334239188212294, 0.18584209157131101)); // 4
     outline.add(new Vec2(3.6497584392382976, 0.18572669434389927)); // 5
     outline.add(new Vec2(4.3086984061600848, 0.18603552950492047)); // 6
     outline.add(new Vec2(5.6832870081551921, 0.18624859438075031)); // 7
     outline.add(new Vec2(5.9833598523943792, 0.18551919930533811)); // 8
     outline.add(new Vec2(6.3177115233380352, 0.18457940020320102)); // 9
     outline.add(new Vec2(6.717001571821406, 0.18524022572244839)); // 10
     outline.add(new Vec2(6.9187933479159405, 0.18707617819954853)); // 11
     outline.add(new Vec2(7.1215523623937713, 0.19036980997708464)); // 12
     outline.add(new Vec2(7.3213353721103545, 0.19534864544021627)); // 13
     outline.add(new Vec2(7.3253214415656345, 0.19546758728278774)); // 14
     outline.add(new Vec2(7.5352189485243564, 0.20290911868504771)); // 15
     outline.add(new Vec2(7.7534552788890156, 0.21326422321622965)); // 16
     outline.add(new Vec2(7.9759153915626495, 0.22684503099235467)); // 17
     outline.add(new Vec2(8.2037333759429298, 0.24429264570141485)); // 18
     outline.add(new Vec2(8.4363411613055206, 0.26640749803671282)); // 19
     outline.add(new Vec2(8.5882065572673341, 0.28354718952780189)); // 20
     outline.add(new Vec2(8.6680596159090726, 0.29351316030470342)); // 21
     outline.add(new Vec2(8.8960056478699681, 0.32598626100198191)); // 22
     outline.add(new Vec2(9.118153827115469, 0.36412887698863228)); // 23
     outline.add(new Vec2(9.3330602130789728, 0.40817289499987724)); // 24
     outline.add(new Vec2(9.5398054084891797, 0.45830715556243257)); // 25
     outline.add(new Vec2(9.7108194374944663, 0.50646526366781097)); // 26
     outline.add(new Vec2(9.7378014954758996, 0.51467827591982873)); // 27
     outline.add(new Vec2(9.9327359296297963, 0.57944001938363288)); // 28
     outline.add(new Vec2(10.129193051457904, 0.65511806790755134)); // 29
     outline.add(new Vec2(10.3219534482391, 0.74045952868119591)); // 30
     outline.add(new Vec2(10.504989185951249, 0.83207281441556558)); // 31
     outline.add(new Vec2(10.67849024007298, 0.92869808824554223)); // 32
     outline.add(new Vec2(10.748222335552903, 0.97028292927689463)); // 33
     outline.add(new Vec2(10.842300413985953, 1.0289405096889359)); // 34
     outline.add(new Vec2(10.998740148909542, 1.1330584513957589)); // 35
     outline.add(new Vec2(11.107348885631321, 1.2102323108904778)); // 36
     outline.add(new Vec2(11.78765958536165, 1.6730529620258834)); // 37
     outline.add(new Vec2(12.45164835000236, 2.1589999999999998)); // 38
     outline.add(new Vec2(-0.1885323943594, 2.1589999999999994)); // 39
  
  theHull.add(new BuoyantBody(300, 300, outline, beam));

  // suss out polygon area by using a MassData object
  // by using a density of 1.0 - the computed "mass" will equal the polygon area.
  mData = new MassData();  
  ((PolygonShape)theHull.get(theHull.size()-1).body.getFixtureList().getShape()).computeMass(mData, 1.0); 
  polyArea = mData.mass;
  
  // Slice mass is its fractional contribution of total hull volume
  theHull.get(theHull.size()-1).setBuoyantBodyMass( ((polyArea*beam)/targetVolume) * targetMass ); // CG will be reset to center of area by default.
  theHull.get(theHull.size()-1).setBuoyantBodyCG(boatCG); // be sure to set CG after setting the mass



// Slice #4
  outline.clear();
  
     outline.add(new Vec2(-0.021638009160505532, 0.24779010485905734)); // 0
     outline.add(new Vec2(0.059748635858484632, 0.24780507581271111)); // 1
     outline.add(new Vec2(1.0800959973470166, 0.24780865107756272)); // 2
     outline.add(new Vec2(1.728730857745385, 0.24784375239517636)); // 3
     outline.add(new Vec2(2.3265217848763218, 0.24778945542841474)); // 4
     outline.add(new Vec2(3.3222386956588257, 0.2477171082150558)); // 5
     outline.add(new Vec2(3.9678660478854697, 0.2480473726732273)); // 6
     outline.add(new Vec2(5.3145025978702538, 0.24808852660206096)); // 7
     outline.add(new Vec2(5.6083929170360909, 0.24733473933376054)); // 8
     outline.add(new Vec2(5.9358482256654295, 0.24641865035976243)); // 9
     outline.add(new Vec2(6.3499012367839693, 0.24710549287980721)); // 10
     outline.add(new Vec2(6.5597777107990929, 0.24901514636621094)); // 11
     outline.add(new Vec2(6.7709263553538364, 0.25244523431045207)); // 12
     outline.add(new Vec2(6.9792401522393277, 0.25763688692933717)); // 13
     outline.add(new Vec2(6.9833973453912304, 0.25776093525757438)); // 14
     outline.add(new Vec2(7.2040082382304309, 0.2655849508734997)); // 15
     outline.add(new Vec2(7.4356993012255934, 0.27658026344278425)); // 16
     outline.add(new Vec2(7.6729768861094918, 0.29106627002819119)); // 17
     outline.add(new Vec2(7.9173657580577599, 0.30978594434565232)); // 18
     outline.add(new Vec2(8.1681994265501014, 0.33363447701367033)); // 19
     outline.add(new Vec2(8.3320122344371654, 0.35212250272901302)); // 20
     outline.add(new Vec2(8.4180162258269, 0.36285619269215719)); // 21
     outline.add(new Vec2(8.6630825432721679, 0.39776615762787992)); // 22
     outline.add(new Vec2(8.900808315351096, 0.43858003588558842)); // 23
     outline.add(new Vec2(9.1293787978703893, 0.48541990469688606)); // 24
     outline.add(new Vec2(9.3476786376971592, 0.53835044634869966)); // 25
     outline.add(new Vec2(9.5268092683860743, 0.58878880896315611)); // 26
     outline.add(new Vec2(9.5550344542333256, 0.59738004541704148)); // 27
     outline.add(new Vec2(9.7591400314952779, 0.66519328123530741)); // 28
     outline.add(new Vec2(9.9662183333879071, 0.74497062465974506)); // 29
     outline.add(new Vec2(10.169404934184339, 0.83491952942571679)); // 30
     outline.add(new Vec2(10.360694505112807, 0.93065618688874019)); // 31
     outline.add(new Vec2(10.540356544741687, 1.0307042048674162)); // 32
     outline.add(new Vec2(10.611807081308395, 1.0733109110843353)); // 33
     outline.add(new Vec2(10.708198988843385, 1.1334071863575859)); // 34
     outline.add(new Vec2(10.733686325965156, 1.1498055318239175)); // 35
     outline.add(new Vec2(12.246206641779999, 2.1590000000000003)); // 36
     outline.add(new Vec2(-0.1885323943594, 2.1589999999999998)); // 37
  
  theHull.add(new BuoyantBody(300, 300, outline, beam));

  // suss out polygon area by using a MassData object
  // by using a density of 1.0 - the computed "mass" will equal the polygon area.
  mData = new MassData();  
  ((PolygonShape)theHull.get(theHull.size()-1).body.getFixtureList().getShape()).computeMass(mData, 1.0); 
  polyArea = mData.mass;
  
  // Slice mass is its fractional contribution of total hull volume
  theHull.get(theHull.size()-1).setBuoyantBodyMass( ((polyArea*beam)/targetVolume) * targetMass ); // CG will be reset to center of area by default.
  theHull.get(theHull.size()-1).setBuoyantBodyCG(boatCG); // be sure to set CG after setting the mass
  

// Slice #5
  outline.clear();
  
     outline.add(new Vec2(-0.027047511450631914, 0.3097376310738218)); // 0
     outline.add(new Vec2(0.074685794823105708, 0.30975634476588892)); // 1
     outline.add(new Vec2(0.90585370660177622, 0.30976081384695342)); // 2
     outline.add(new Vec2(1.4342259595392288, 0.30978885491121333)); // 3
     outline.add(new Vec2(2.0196196509314124, 0.30973681928551855)); // 4
     outline.add(new Vec2(2.9947189520793533, 0.30970752208621199)); // 5
     outline.add(new Vec2(3.6270336896108537, 0.31005921584153401)); // 6
     outline.add(new Vec2(4.9457181875853156, 0.30992845882337161)); // 7
     outline.add(new Vec2(5.2334259816778044, 0.30915027936218276)); // 8
     outline.add(new Vec2(5.5539849279928237, 0.3082579005163239)); // 9
     outline.add(new Vec2(5.9828009017465327, 0.30897076003716606)); // 10
     outline.add(new Vec2(6.2007620736822471, 0.31095411453287325)); // 11
     outline.add(new Vec2(6.4203003483139032, 0.31452065864381962)); // 12
     outline.add(new Vec2(6.6371449323683009, 0.31992512841845822)); // 13
     outline.add(new Vec2(6.6414732492168254, 0.32005428323236107)); // 14
     outline.add(new Vec2(6.8727975279365054, 0.32826078306195183)); // 15
     outline.add(new Vec2(7.1179433235621739, 0.33989630366933876)); // 16
     outline.add(new Vec2(7.3700383806563341, 0.35528750906402795)); // 17
     outline.add(new Vec2(7.6309981401725864, 0.37527924298988963)); // 18
     outline.add(new Vec2(7.9000576917946841, 0.4008614559906275)); // 19
     outline.add(new Vec2(8.0758179116070004, 0.42069781593022415)); // 20
     outline.add(new Vec2(8.1679728357447292, 0.43219922507961117)); // 21
     outline.add(new Vec2(8.4301594386743659, 0.46954605425377777)); // 22
     outline.add(new Vec2(8.6834628035867212, 0.51303119478254455)); // 23
     outline.add(new Vec2(8.9256973826618058, 0.56266691439389516)); // 24
     outline.add(new Vec2(9.1555518669051388, 0.61839373713496659)); // 25
     outline.add(new Vec2(9.3427990992776806, 0.67111235425850135)); // 26
     outline.add(new Vec2(9.3722674129907499, 0.68008181491425479)); // 27
     outline.add(new Vec2(9.5855441333607594, 0.75094654308698194)); // 28
     outline.add(new Vec2(9.8032436153179088, 0.8348231814119389)); // 29
     outline.add(new Vec2(10.016856420129576, 0.92937953017023744)); // 30
     outline.add(new Vec2(10.216399824274363, 1.0292395593619146)); // 31
     outline.add(new Vec2(10.317437423495381, 1.0842151006319101)); // 32
     outline.add(new Vec2(10.777279428342036, 1.3813410507230257)); // 33
     outline.add(new Vec2(11.992474031214813, 2.1589999999999998)); // 34
     outline.add(new Vec2(-0.1885323943594, 2.1589999999999998)); // 35
  
  theHull.add(new BuoyantBody(300, 300, outline, beam));

  // suss out polygon area by using a MassData object
  // by using a density of 1.0 - the computed "mass" will equal the polygon area.
  mData = new MassData();  
  ((PolygonShape)theHull.get(theHull.size()-1).body.getFixtureList().getShape()).computeMass(mData, 1.0); 
  polyArea = mData.mass;
  
  // Slice mass is its fractional contribution of total hull volume
  theHull.get(theHull.size()-1).setBuoyantBodyMass( ((polyArea*beam)/targetVolume) * targetMass ); // CG will be reset to center of area by default.
  theHull.get(theHull.size()-1).setBuoyantBodyCG(boatCG); // be sure to set CG after setting the mass


// Slice #6
  outline.clear();
  
     outline.add(new Vec2(-0.032457013740758406, 0.37168515728858609)); // 0
     outline.add(new Vec2(0.089622953787726778, 0.37170761371906669)); // 1
     outline.add(new Vec2(0.731611415856536, 0.37171297661634423)); // 2
     outline.add(new Vec2(1.139721061333073, 0.37173395742725046)); // 3
     outline.add(new Vec2(1.712717516986503, 0.37168418314262219)); // 4
     outline.add(new Vec2(2.6671992084998806, 0.37169793595736855)); // 5
     outline.add(new Vec2(3.2862013313362377, 0.37207105900984105)); // 6
     outline.add(new Vec2(4.5769337773003773, 0.37176839104468218)); // 7
     outline.add(new Vec2(4.8584590463195161, 0.37096581939060508)); // 8
     outline.add(new Vec2(5.1721216303202189, 0.37009715067288518)); // 9
     outline.add(new Vec2(5.615700566709096, 0.37083602719452474)); // 10
     outline.add(new Vec2(5.8417464365654013, 0.37289308269953603)); // 11
     outline.add(new Vec2(6.06967434127397, 0.37659608297718727)); // 12
     outline.add(new Vec2(6.2950497124972724, 0.38221336990757948)); // 13
     outline.add(new Vec2(6.299549153042423, 0.38234763120714765)); // 14
     outline.add(new Vec2(6.5415868176425791, 0.39093661525040413)); // 15
     outline.add(new Vec2(6.8001873458987507, 0.40321234389589344)); // 16
     outline.add(new Vec2(7.0670998752031764, 0.41950874809986455)); // 17
     outline.add(new Vec2(7.3446305222874138, 0.44077254163412705)); // 18
     outline.add(new Vec2(7.6319159570392667, 0.4680884349675849)); // 19
     outline.add(new Vec2(7.8196235887768335, 0.48927312913143534)); // 20
     outline.add(new Vec2(7.9179294456625549, 0.50154225746706493)); // 21
     outline.add(new Vec2(8.1972363340765639, 0.54132595087967572)); // 22
     outline.add(new Vec2(8.4661172918223464, 0.58748235367950086)); // 23
     outline.add(new Vec2(8.7220159674532241, 0.63991392409090408)); // 24
     outline.add(new Vec2(8.9634250961131166, 0.69843702792123352)); // 25
     outline.add(new Vec2(9.1587889301692869, 0.75343589955384638)); // 26
     outline.add(new Vec2(9.1895003717481742, 0.76278358441146765)); // 27
     outline.add(new Vec2(9.4119482352262409, 0.83669980493865659)); // 28
     outline.add(new Vec2(9.6402688972479105, 0.92467573816413262)); // 29
     outline.add(new Vec2(9.837122839276276, 1.011104131419359)); // 30
     outline.add(new Vec2(11.19811940528184, 1.8627831322545552)); // 31
     outline.add(new Vec2(11.681288056051642, 2.1590000000000003)); // 32
     outline.add(new Vec2(-0.1885323943594, 2.1589999999999989)); // 33
  
  theHull.add(new BuoyantBody(300, 300, outline, beam));

  // suss out polygon area by using a MassData object
  // by using a density of 1.0 - the computed "mass" will equal the polygon area.
  mData = new MassData();  
  ((PolygonShape)theHull.get(theHull.size()-1).body.getFixtureList().getShape()).computeMass(mData, 1.0); 
  polyArea = mData.mass;
  
  // Slice mass is its fractional contribution of total hull volume
  theHull.get(theHull.size()-1).setBuoyantBodyMass( ((polyArea*beam)/targetVolume) * targetMass ); // CG will be reset to center of area by default.
  theHull.get(theHull.size()-1).setBuoyantBodyCG(boatCG); // be sure to set CG after setting the mass


// Slice #7
  outline.clear();
  
     outline.add(new Vec2(-0.037866516030884687, 0.43363268350335055)); // 0
     outline.add(new Vec2(0.10456011275234797, 0.43365888267224428)); // 1
     outline.add(new Vec2(0.55736912511129566, 0.43366513938573481)); // 2
     outline.add(new Vec2(0.84521616312691727, 0.43367905994328732)); // 3
     outline.add(new Vec2(1.4058153830415938, 0.43363154699972584)); // 4
     outline.add(new Vec2(2.3396794649204091, 0.43368834982852483)); // 5
     outline.add(new Vec2(2.9453689730616213, 0.43408290217814771)); // 6
     outline.add(new Vec2(4.2081493670154391, 0.43360832326599286)); // 7
     outline.add(new Vec2(4.4834921109612278, 0.43278135941902729)); // 8
     outline.add(new Vec2(4.7902583326476122, 0.43193640082944673)); // 9
     outline.add(new Vec2(5.2486002316716593, 0.43270129435188343)); // 10
     outline.add(new Vec2(5.4827307994485537, 0.43483205086619814)); // 11
     outline.add(new Vec2(5.7190483342340359, 0.4386715073105546)); // 12
     outline.add(new Vec2(5.9529544926262457, 0.44450161139670019)); // 13
     outline.add(new Vec2(5.957625056868018, 0.44464097918193435)); // 14
     outline.add(new Vec2(6.2103761073486545, 0.45361244743885626)); // 15
     outline.add(new Vec2(6.4824313682353303, 0.4665283841224479)); // 16
     outline.add(new Vec2(6.7641613697500187, 0.48372998713570098)); // 17
     outline.add(new Vec2(7.0582629044022429, 0.50626584027836408)); // 18
     outline.add(new Vec2(7.3637742222838476, 0.5353154139445423)); // 19
     outline.add(new Vec2(7.5634292659466675, 0.55784844233264641)); // 20
     outline.add(new Vec2(7.6678860555803841, 0.57088528985451914)); // 21
     outline.add(new Vec2(7.9643132294787637, 0.61310584750557351)); // 22
     outline.add(new Vec2(8.2487717800579716, 0.66193351257645716)); // 23
     outline.add(new Vec2(8.5183345522446423, 0.71716093378791301)); // 24
     outline.add(new Vec2(8.7712983253210961, 0.77848031870750045)); // 25
     outline.add(new Vec2(8.9747787610608949, 0.83575944484919162)); // 26
     outline.add(new Vec2(9.0067333305055985, 0.8454853539086804)); // 27
     outline.add(new Vec2(9.2383523370917207, 0.92245306679033146)); // 28
     outline.add(new Vec2(9.249076540626131, 0.92630207234910844)); // 29
     outline.add(new Vec2(10.269807313683629, 1.5528727420267845)); // 30
     outline.add(new Vec2(11.302808699117442, 2.1590000000000007)); // 31
     outline.add(new Vec2(-0.1885323943594, 2.1589999999999998)); // 32
  
  theHull.add(new BuoyantBody(300, 300, outline, beam));

  // suss out polygon area by using a MassData object
  // by using a density of 1.0 - the computed "mass" will equal the polygon area.
  mData = new MassData();  
  ((PolygonShape)theHull.get(theHull.size()-1).body.getFixtureList().getShape()).computeMass(mData, 1.0); 
  polyArea = mData.mass;
  
  // Slice mass is its fractional contribution of total hull volume
  theHull.get(theHull.size()-1).setBuoyantBodyMass( ((polyArea*beam)/targetVolume) * targetMass ); // CG will be reset to center of area by default.
  theHull.get(theHull.size()-1).setBuoyantBodyCG(boatCG); // be sure to set CG after setting the mass


// Slice #8
  outline.clear();
  
     outline.add(new Vec2(-0.043276018321011064, 0.49558020971811501)); // 0
     outline.add(new Vec2(0.11949727171696904, 0.49561015162542232)); // 1
     outline.add(new Vec2(0.38312683436605544, 0.49561730215512556)); // 2
     outline.add(new Vec2(0.55071126492076128, 0.49562416245932434)); // 3
     outline.add(new Vec2(1.0989132490966851, 0.49557891085682959)); // 4
     outline.add(new Vec2(2.0121597213409368, 0.49567876369968106)); // 5
     outline.add(new Vec2(2.6045366147870044, 0.4960947453464547)); // 6
     outline.add(new Vec2(3.8393649567305017, 0.49544825548730342)); // 7
     outline.add(new Vec2(4.1085251756029404, 0.49459689944744956)); // 8
     outline.add(new Vec2(4.4083950349750083, 0.49377565098600795)); // 9
     outline.add(new Vec2(4.8814998966342209, 0.49456656150924239)); // 10
     outline.add(new Vec2(5.1237151623317079, 0.49677101903286069)); // 11
     outline.add(new Vec2(5.3684223271941018, 0.50074693164392237)); // 12
     outline.add(new Vec2(5.6108592727552189, 0.50678985288582146)); // 13
     outline.add(new Vec2(5.6157009606936157, 0.50693432715672104)); // 14
     outline.add(new Vec2(5.8791653970547273, 0.51628827962730828)); // 15
     outline.add(new Vec2(6.1646753905719081, 0.52984442434900236)); // 16
     outline.add(new Vec2(6.4612228642968592, 0.5479512261715378)); // 17
     outline.add(new Vec2(6.7718952865170712, 0.57175913892260144)); // 18
     outline.add(new Vec2(7.0956324875284293, 0.60254239292149947)); // 19
     outline.add(new Vec2(7.3072349431164998, 0.62642375553385765)); // 20
     outline.add(new Vec2(7.4178426654982097, 0.64022832224197312)); // 21
     outline.add(new Vec2(7.7313901248809636, 0.68488574413147152)); // 22
     outline.add(new Vec2(8.0314262682936004, 0.73638467147341347)); // 23
     outline.add(new Vec2(8.314653137036057, 0.79440794348492205)); // 24
     outline.add(new Vec2(8.4101074543444501, 0.81634137516239946)); // 25
     outline.add(new Vec2(9.0552776566756386, 1.1847639014335005)); // 26
     outline.add(new Vec2(10.834465291334624, 2.1590000000000003)); // 27
     outline.add(new Vec2(-0.1885323943594, 2.1589999999999994)); // 28
  
  theHull.add(new BuoyantBody(300, 300, outline, beam));

  // suss out polygon area by using a MassData object
  // by using a density of 1.0 - the computed "mass" will equal the polygon area.
  mData = new MassData();  
  ((PolygonShape)theHull.get(theHull.size()-1).body.getFixtureList().getShape()).computeMass(mData, 1.0); 
  polyArea = mData.mass;
  
  // Slice mass is its fractional contribution of total hull volume
  theHull.get(theHull.size()-1).setBuoyantBodyMass( ((polyArea*beam)/targetVolume) * targetMass ); // CG will be reset to center of area by default.
  theHull.get(theHull.size()-1).setBuoyantBodyCG(boatCG); // be sure to set CG after setting the mass




// Slice #9
  outline.clear();
  
//#9
     outline.add(new Vec2(-0.096362146030628301, 1.1035019948310454)); // 0
     outline.add(new Vec2(1.7521705574102993, 0.9801361131633981)); // 1
     outline.add(new Vec2(3.5805466324053032, 0.89991419035939368)); // 2
     outline.add(new Vec2(5.3980729982183391, 0.85499514610131111)); // 3
     outline.add(new Vec2(5.8609624904570889, 0.83827317219658148)); // 4
     outline.add(new Vec2(6.324088978709697, 0.84602960920889192)); // 5
     outline.add(new Vec2(6.5580997711355682, 0.86731957224142053)); // 6
     outline.add(new Vec2(6.7844968948192932, 0.90085577474572975)); // 7
     outline.add(new Vec2(7.0088398568570849, 0.94630819602212957)); // 8
     outline.add(new Vec2(7.2366881643451615, 1.0033468153709295)); // 9
     outline.add(new Vec2(7.7869760220477922, 1.1726911086320229)); // 10
     outline.add(new Vec2(8.3292462622479704, 1.3661780657800588)); // 11
     outline.add(new Vec2(9.3968356772997694, 1.7969710026253831)); // 12
     outline.add(new Vec2(10.194551179648787, 2.1589999999999998)); // 13
     outline.add(new Vec2(-0.1885323943594, 2.1589999999999994)); // 14
  
  theHull.add(new BuoyantBody(305, 300, outline, beam));

  // suss out polygon area by using a MassData object
  // by using a density of 1.0 - the computed "mass" will equal the polygon area.
  mData = new MassData();  
  ((PolygonShape)theHull.get(theHull.size()-1).body.getFixtureList().getShape()).computeMass(mData, 1.0); 
  polyArea = mData.mass;
  
  // Slice mass is its fractional contribution of total hull volume
  theHull.get(theHull.size()-1).setBuoyantBodyMass( ((polyArea*beam)/targetVolume) * targetMass ); // CG will be reset to center of area by default.
  theHull.get(theHull.size()-1).setBuoyantBodyCG(boatCG); // be sure to set CG after setting the mass



  
  // Create the new HullBody.
  hull = new HullBody(theHull);
  
  

  // Create PID Controller
  PID = new PIDController();

  PID.iT = PID.dT = dt;  // set time constants to the world time step;
  PID.kP = 1.0;          // Proportional controller constant
  PID.kI = 1.0;          // Integral controller constant
  PID.kD = 0.0;          // Derivative controller constant
  
  PID.setpoint=15.0;  // target speed in m/s


  // Set up data file output.
  output = createWriter("accels.txt"); 
  
}



void draw() {

  // Advance the wave
  sigWaveHt=sldrSigWaveHt.getValueF();
  surface.SetWaveParams(6, sigWaveHt);  // Wave period (s), H13 (m)
  surface.AdvanceWave(t);

  // Use PID controler to regulate speed
  PID.setpoint=sldrTargetVelocity.getValueF();  // Get target velocity from GUI control
  float throttleSetting = PID.PIDOutputValue(hull.buoyantBodies.get(0).body.getLinearVelocity().x);  // Compute throttle setting from VMG
  // Apply thrust to first hull body (TODO: add shaft angle for thrust vector)
  hull.buoyantBodies.get(0).setBuoyantBodyThrustVector(new Vec2(maxThrust,0.0).mul(throttleSetting)); 

  // step time through the world
  box2d.step();

  scale(1);
  //translate(-500/2,-300/2);

  // Translate the viewpoint to travel with the boat
  pushMatrix();
  
  Vec2 v = box2d.coordWorldToPixels(hull.buoyantBodies.get(0).body.getWorldCenter());
  translate((v.x * -1)+(width/2),0);

// Render the scene
  background(255);
  
  // Draw the tank
  tank.display();
  
  // Draw the surface
  surface.display();
  
  // Draw the boat;
  //boat.display();
  //boat.processBuoyancy();

  // Draw the hull;
  hull.processAndDisplay();
  
  popMatrix();

  Vec2 currentVel = new Vec2(hull.buoyantBodies.get(0).body.getLinearVelocity());
  dv.set(prevV.subLocal(currentVel));
  prevV.set(currentVel);
  
  Vec2 gravity = box2d.world.getGravity();
  Vec2 accel= new Vec2(dv).mulLocal(1.0 / dt);
  
// and time marches on...
  t += dt;



// ********************** On-screen Stats ************************ //

  text("Trimming Angle (deg):" + hull.buoyantBodies.get(0).body.getAngle()*(180/PI), 300,128);
  text("Speed (VMG) (m/s):" + currentVel.x + "  Knots:" + currentVel.x*1.9438 , 300,144);
  text("Velocity (Vector) (m/s):" + currentVel, 300,160);
  text("dVelocity (m/s):" + dv, 300,176);
  text("Acceleration (m/s^2):" + accel, 300,192);
  text("Acceleration (g): " + accel.x / gravity.length() + "," + accel.y / gravity.length(), 300,208);
  text("BoatBodyMass (kg):" + hull.buoyantBodies.get(0).body.getMass(), 300,224);
  text("Throttle Setting (%):" + throttleSetting, 300,240);
  text("Velocity Target (m/s):" + PID.setpoint, 300,250);
  
  
  // Send data to outFile
  output.println(
                t + "," + 
                throttleSetting+ "," +
                currentVel.x+ "," + 
                hull.buoyantBodies.get(0).body.getAngle()*(180/PI) + ","  + 
                accel.x / gravity.length() + "," + 
                accel.y / gravity.length() 
                );


  // Framerate and other stats
  fill(0);
  text("framerate: " + (int)frameRate,12,16);
  
  Vec2 mousePos= new Vec2(mouseX, mouseY);
  text("Mouse (window):" + mousePos,300,16);
  text("Mouse (world) :" + box2d.coordPixelsToWorld(mousePos) ,300,32);
  
  text("Elasped Simulation Time (s):" + t ,500,16);
  text("Wave Period (s):" + surface.wavePeriod ,500,32);
  text("Sig. Wave Height [h13] (m):" + surface.h13 ,500,48);
  text("Wavelength (m):" + surface.waveLength ,500,64);
  text("Wave speed (m/s):" + surface.waveSpeed ,500,80);

  
//  delay(125);
//  saveFrame("frames/####.png");
}



// Collision event functions!
void beginContact(Contact cp) {
  // Get both fixtures
  Fixture f1 = cp.getFixtureA();
  Fixture f2 = cp.getFixtureB();

  // Get both bodies from thier fixtures
  Body b1 = f1.getBody();
  Body b2 = f2.getBody();
  
  // Get the objects that reference these bodies from thier fixtures (requires "body.setUserData(this);" during body creation) 
  Object o1 = b1.getUserData();
  Object o2 = b2.getUserData();
  
  BuoyantBody bb = null;
  WaveBody wb = null;

  // We don't know the order or class of contacting bodies (yet), so lets sort that out.
  if (o1.getClass() == BuoyantBody.class && o2.getClass() == WaveBody.class) {
    bb = (BuoyantBody)o1;
    wb = (WaveBody)o2;
    
  }
  else if (o1.getClass() == WaveBody.class && o2.getClass() == BuoyantBody.class) {
    bb = (BuoyantBody)o2;
    wb = (WaveBody)o1;
  }

  // Now tell the BuoyantBody which WaveBody it's associated with.
  if(bb != null && wb != null){
    bb.change();
    bb.setContactWaveBodyObject(wb);
  }
}





// Objects stop touching each other
void endContact(Contact cp) {
  // Get both fixtures
  Fixture f1 = cp.getFixtureA();
  Fixture f2 = cp.getFixtureB();

  // Get both bodies from thier fixtures
  Body b1 = f1.getBody();
  Body b2 = f2.getBody();

  // Get the objects that reference these bodies
  Object o1 = b1.getUserData();
  Object o2 = b2.getUserData();
  
  BuoyantBody bb = null;
  WaveBody wb = null;

  // We don't know the order or class of contacting bodies (yet), so lets sort that out.
  if (o1.getClass() == BuoyantBody.class && o2.getClass() == WaveBody.class) {
    bb = (BuoyantBody)o1;
    wb = (WaveBody)o2;
  }
else if(o1.getClass() == WaveBody.class && o2.getClass() == BuoyantBody.class){
    bb = (BuoyantBody)o2;
    wb = (WaveBody)o1;
  }
  
  // Now tell the BuoyantBody to clear the WaveBody it's currently associated with.
  // If the Buoyant Body is contained entirely within the Wavebody, the association will not be cleared.
    if(bb != null && wb != null){
      bb.clearContactWaveBodyObject();
    }
    
}



/* 
void mousePressed() {
boolean toggle = false;

if(toggle){
  toggle=false;
  loop();
  }

else{
  toggle=true;
  noLoop();
}

}
*/


/*
void mouseReleased() {
  loop();
}
*/

public void handleSliderEvents(GValueControl slider, GEvent event) { 
//  if (slider == sdr)  // The slider being configured?
//    println(sdr.getValueS() + "    " + event);    

  if (slider == sldrSigWaveHt){
    //sdr.setEasing(slider.getValueF());
    sigWaveHt=sldrSigWaveHt.getValueF();
  }

  else if (slider == sldrTargetVelocity){
    //sdr.setEasing(slider.getValueF());
    PID.setpoint=sldrTargetVelocity.getValueF();  // target speed in m/s
  }

}