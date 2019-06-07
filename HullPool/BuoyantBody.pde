class BuoyantBody {

  // We need to keep track of a Body and a shape
  Body body; 
  PolygonShape cs;

  //Body coordinates and color
  ArrayList<Vec2> coords;
  color col;

  // Mass of the shape (kg)
  float kg;

/* beam of the fixture is in FixtureDef.userData  
  // Beam (width) of the shape (m);
  float b;
*/

  // Center of gravity - relative to shape origin
  Vec2 cg;
  
  // Thrust
  Vec2 propulsorPoint; // Location of the propulsor - relative to the shape origin
  Vec2 thrustVector; //  Vector of propulsion thrust - relative to the shape (origin)
  

  // The WaveBody being contacted.  null if no contact
  WaveBody contactWaveBodyObject=null;


// New Buoyant Body
// x & y are in screen (pixel) coordinates
// coords_ are in world coordinates
  BuoyantBody(float x, float y, ArrayList<Vec2> coords_, float beam) {
    coords = new ArrayList<Vec2>(coords_);


    // This function puts the particle in the Box2d world
    makeBody(x, y, coords, beam);
    body.setUserData(this);
    
    setBuoyantBodyCG(body.getLocalCenter());  // Default CG is center of Area - relative to shape origin
    
    setBuoyantBodyPropulsorPoint(body.getLocalCenter());  // Default thust point is through the center of area (default CG);
    setBuoyantBodyThrustVector(new Vec2(0,0));  // Zero thrust is the default

    setBuoyantBodyMass(1f);  // Default mass (mass is used in makeBody()
 
    

    col = color(175);
  }



// Set the buoyant body mass (kg)
// and fix up body density.
void setBuoyantBodyMass(float m_){
  kg = m_;
  
  Shape shape;
  Fixture fixture;
  
  fixture=body.getFixtureList();  //get the fixture
  
  shape=fixture.getShape(); // get shape associated with fixture

  // suss out polygon area by using a MassData object
  MassData mData = new MassData();
  cs.computeMass(mData, 1.0); // by using a density of 1.0 - the computed "mass" will equal the polygon area.
  float polyArea = mData.mass;
      
  // set density of the fixture.
  fixture.setDensity(kg/polyArea);

  
  body.resetMassData(); //update the body mass data
  
}





/*
// Set the buoyant body beam (m)
void setBuoyantBodyBeam(float b_){
  b = b_;
}
*/

// Set the propulsion thrust point
void setBuoyantBodyPropulsorPoint(Vec2 tPoint){
  propulsorPoint = new Vec2(tPoint);
}

// Set the propulsion thrust vector
void setBuoyantBodyThrustVector(Vec2 tVec){
  thrustVector = new Vec2(tVec);
}


// Set the true Center of Gravity
void setBuoyantBodyCG(Vec2 _cg){
  cg = new Vec2(_cg);
  
  MassData data = new MassData();
  body.getMassData(data); //get the mass data from the body

  /** Set the position of the shape's centroid relative to the shape's origin. **/
  data.center.set(cg.x, cg.y);

  body.setMassData(data); //update the body's mass data
}


// This function removes the particle from the box2d world
  void killBody() {
    box2d.destroyBody(body);
  }


// **
// Change color when hit
//
  void change() {
    col = color(255, 0, 0);
  }


// Is the particle ready for deletion?
  boolean done() {
    // Let's find the screen position of the particle
    Vec2 pos = box2d.getBodyPixelCoord(body);
    // Is it off the bottom of the screen?
    if (pos.y > height /* +r*2 */) {
      killBody();
      return true;
    }
    return false;
  }



// Here's our function that adds the body to the Box2D world
  void makeBody(float x, float y, ArrayList<Vec2> coords, float beam) {
    // Define the body
    BodyDef bd = new BodyDef();
    // Set its position
    bd.position = box2d.coordPixelsToWorld(x, y);
    bd.type = BodyType.DYNAMIC;
//  bd.bullet=true;

    // Create Body
    body = box2d.createBody(bd);

    // Define the Shape
    // Make the body's shape a polygon    
    cs = new PolygonShape();

    Vec2[] vertices = new Vec2[coords.size()];
    for(int i=0; i < coords.size(); ++i)
      vertices[i]=coords.get(i);
        
    cs.set(vertices, vertices.length);

    // compute shape centroid in world coordinates
    Vec2 centroid = new Vec2();
    cs.computeCentroidToOut(vertices, vertices.length, centroid);
    
    // suss out polygon area by using a MassData object
    MassData mData = new MassData();
    cs.computeMass(mData, 1.0); // by using a density of 1.0 - the computed "mass" will equal the polygon area.
    float polyArea = mData.mass;
    

    




    // Define the Fixture
    FixtureDef fd = new FixtureDef();

    fd.filter.groupIndex = -1;  // use a negative group index to turn off collision processing between BuouantBody fixtures
                                // Note:  This worked where 'joint.collideConnected = false' did not. 

    fd.shape = cs;
    // Parameters that affect physics
    fd.density = kg/polyArea;
    
    fd.friction = 0.01;
    fd.restitution = 0.3;
    
    fd.userData = beam;

    // Attach fixture to body
    body.createFixture(fd);
    body.setAngularVelocity(0);
    

  }


// Set our WaveBody when we contact it
// Called from beginContact callback function
  void setContactWaveBodyObject(WaveBody wbo){
    contactWaveBodyObject = wbo;
  }


// Clear the WaveBody when it's no longer in contact with BuoyantBody - or Buoyant Body is contained entirely within Wavebody
// Called from endContact callback function
 void clearContactWaveBodyObject(){
 
  if(contactWaveBodyObject != null){
  
    // Water geometry - in its local world frame
    Body waveBody = contactWaveBodyObject.body;
    Fixture waveFixture = waveBody.getFixtureList();
     
    Vec2 waterWorldVertices[] = new Vec2[((ChainShape)waveFixture.getShape()).m_count];
    for(int i=0; i < waterWorldVertices.length; ++i)
      waterWorldVertices[i] = waveBody.getWorldPoint(((ChainShape)waveFixture.getShape()).m_vertices[i]);    
    
    
    boolean inWater=true;
    for(int i=0; i < coords.size(); ++i){
      Vec2 bPointWorld = body.getWorldPoint(coords.get(i)); // body point in current world coordinates
      
      if(pointInPoly(bPointWorld, waterWorldVertices, waterWorldVertices.length) != 0) { // otherwise test to see if body point is in the water
        inWater=true;
      }
      else {
        inWater=false;
        break;
      }
    }
    
    if(inWater==false){  
      contactWaveBodyObject = null;
    }
  
  }

}


// isLeft and pointInPoly adapted from:
// http://geomalgorithms.com/a03-_inclusion.html

// isLeft: tests if a point is Left|On|Right of an infinite line.
//    Input:  three points A, B, and pt
//    Return: >0 for pt left of the line through A and B
//            =0 for pt  on the line
//            <0 for pt  right of the line
private float isLeft(Vec2 a, Vec2 b, Vec2 pt )
{
    return ( (b.x - a.x) * (pt.y - a.y) - (pt.x -  a.x) * (b.y - a.y) );
}


// pointInPoly: winding number test for a point in a polygon
//      Input:   P = the point to test,
//               V[] = vertex points of the polygon
//      Return:  wn = the winding number (Equal to 0 when P is outside)
int pointInPoly(Vec2 P, Vec2[] V, int n )
{
    int    wn = 0;    // the  winding number counter
 
    // loop through all edges of the polygon
    for (int i=0; i < n; i++) {   // edge from V[i] to  V[i+1]
        int nextIndex = (i+1) % n;
        if (V[i].y <= P.y) {          // start y <= P.y
            if (V[nextIndex].y  > P.y)      // an upward crossing
                 if (isLeft( V[i], V[nextIndex], P) > 0)  // P left of  edge
                     ++wn;            // have  a valid up intersect
        }
        else {                        // start y > P.y (no test needed)
            if (V[nextIndex].y  <= P.y)     // a downward crossing
                 if (isLeft( V[i], V[nextIndex], P) < 0)  // P right of  edge
                     --wn;            // have  a valid down intersect
        }
    }

    return(wn);
}

/********************** BEGIN BUOYANCY *************************************/

void processBuoyancy(){

  // If we are not in contact with a WaveBody, what are we doing here?  
  if (contactWaveBodyObject == null)
    return;

  Body waveBody = contactWaveBodyObject.body;
  Fixture waveFixture = waveBody.getFixtureList();
     
  // Stuff about the world  
  float density = waveFixture.getDensity();
  Vec2 gravity = box2d.world.getGravity();
      
 
  // Translate local shapes to current world coordinates
  // Buoyant geometry - in its local world frame
  Vec2 buoyWorldVertices[] = new Vec2[coords.size()];
  for(int i=0; i < buoyWorldVertices.length; ++i)
    buoyWorldVertices[i] = body.getWorldPoint(coords.get(i));    


  // Water geometry - in its local world frame
     Vec2 waterWorldVertices[] = new Vec2[((ChainShape)waveFixture.getShape()).m_count];
 
  for(int i=0; i < waterWorldVertices.length; ++i)
    waterWorldVertices[i] = waveBody.getWorldPoint(((ChainShape)waveFixture.getShape()).m_vertices[i]);    


  // Subject and Clipping polygon inputs arrays of Vec2
  // Result is provided as an array of Vec2 or ArrayList of Vec2
  // Result geometry and centroid will be in world unit coordinates.  Area will be in world units.
  SutherlandHodgmanClipper submergedPortion = new SutherlandHodgmanClipper(waterWorldVertices, buoyWorldVertices);


  //********************* HYDROSTATIC FORCES ****************************//
  // Get submerged area and centroid
  float area = submergedPortion.resultCenterArea.area;
  Vec2 buoyCentroid = submergedPortion.resultCenterArea.centroid;
  float beam =(Float)body.getFixtureList().m_userData;  
 
  // Apply buoyancy force to the submerged portion our body

  float displacedMass = density * submergedPortion.resultCenterArea.area * beam; 
  Vec2 buoyancyForce = new Vec2().set(gravity).negate().mulLocal(displacedMass);
  body.applyForce( buoyancyForce, submergedPortion.resultCenterArea.centroid );

/*
  // Apply trimming moment to account for actual CG location from center of body area.
  Vec2 cgWorld = body.getWorldPoint(cg);
  Vec2 bodyCenter = new Vec2(body.getWorldCenter()); 
  Vec2 gz = bodyCenter.subLocal(cgWorld);
  float arm=gz.length();
  body.applyTorque(arm * buoyancyForce.length()); 
*/
  
  //********************* PROPULSION FORCES ****************************//
  // Apply propulsion
  Vec2 tPointWorld = body.getWorldPoint(propulsorPoint); // Thrust point in current world coordinates
  
  // Test that propulson point is "in the water" - that is, the propulsor point is within the local water body perimeter.
  // In the special case where the proplusor thrust point is at the boat CG, then thrust is always applied
  boolean inWater;
  
  if(propulsorPoint.x == cg.x && propulsorPoint.y == cg.y) // prolusion point is at the CG
    inWater=true;    
  else if(pointInPoly(tPointWorld, waterWorldVertices, waterWorldVertices.length) != 0) // otherwise test to see if propulsion point is in the water
    inWater=true;
  else
    inWater=false;
    
  if(inWater==true){  
    text("THRUST" ,500,200);


    float theta = body.getAngle(); // Current trim angle
theta=hull.buoyantBodies.get(0).body.getAngle();
    Vec2 rotateThrust = new Vec2(thrustVector.x*cos(theta),thrustVector.y*sin(theta)); // Thrust vector relative to trim angle
    body.applyForce(rotateThrust, tPointWorld);
  }
  else
      text("NO THRUST" ,500,200);
  
  
  //********************* HYDRODYNAMIC FORCES ****************************//
  // Apply dynamic drag and lift separately for each polygon edge
  
    for(int i = 0; i < submergedPortion.result.size(); ++i) {
      //the end points and mid-point of this edge 
      Vec2 p0 = submergedPortion.result.get( (i+0) % submergedPortion.result.size() );
      Vec2 p1 = submergedPortion.result.get( (i+1) % submergedPortion.result.size() );
      Vec2 midPoint = new Vec2(p0).addLocal(p1).mulLocal(0.5f);

      
      //find relative velocity between object and fluid at edge midpoint
      Vec2 velDir1 = waveBody.getLinearVelocityFromWorldPoint( midPoint );
      Vec2 velDir2 = body.getLinearVelocityFromWorldPoint( midPoint );
      Vec2 velDir= new Vec2(velDir2).subLocal(velDir1);
      float vel = velDir.normalize(); // vel has velocity magnitude, velDir is "hat"
  
      Vec2 edge = new Vec2(p1).subLocal(p0);  
      float edgeLength = edge.normalize(); // edgeLength has length, edge is "hat"
      Vec2 normal = Vec2.cross(-1,edge); //gets perpendicular vector (outward pointing normal - assuming ccw body winding)
  
  
      float liftDot = Vec2.dot(edge, velDir); // edge^ dot velDir^ p    
      float dragDot = Vec2.dot(normal, velDir);
      if ( dragDot < 0  || liftDot < 0)
          continue; // Don't process this edge if the normal points backwards relative to velocity - this is not a leading edge
  
  
    stroke(0,0,255);
  strokeWeight(3);
    Vec2 pv0 = box2d.coordWorldToPixels(p0);
        Vec2 pv1 = box2d.coordWorldToPixels(p1);
   line(pv0.x, pv0.y, pv1.x, pv1.y);
      point(pv0.x, pv0.y);
//      point(pv1.x, pv1.y);
      text((i+0) % submergedPortion.result.size() , pv0.x, pv0.y+16);
      text((i+1) % submergedPortion.result.size() , pv0.x, pv0.y+32);


      float nu = 1.19074e-6;
      float cv = vel/sqrt((gravity.length()*beam)); 
      float tau = acos(liftDot);
      float lk = edgeLength;
      float lambda = (lk/beam);
      float v1x = 1.0 - (0.0120*pow(tau*(180.0/PI), 1.1)) / sqrt(lambda);
      float v1 = vel * ((v1x>=0.0) ? sqrt(v1x): 1.0);


      float rn = (vel * lk) / nu;
      float cf = 0.075 / pow( (log(rn)/log(10.0)) - 2.0, 2.0) + 0.0004;
      float cl0 = pow(tau*(180.0/PI), 1.1)*(0.0120 * pow(lambda, 0.5) + 0.0055 * pow(lambda, 5.0/2.0) / pow(cv, 2.0));
  
      // CdP expression fitred data from: Ergebnisse der Aerodynamischen Versuchsanstalt zu Göttingen IV (Prandl - Betz_ 1932)
      //  as cited by Hoerner in Fluid Dynamic Drag

      float cdP = (1.16333 + 0.812234 * exp(-17.9659 * lambda)) - 0.97;  // Factor at end is to correlate full scale test data speed.
      // Drag Force
      float pressureDrag = (0.5 * density * vel * vel) * (lk * beam) * (cdP * dragDot);
      float viscousDrag = (0.5 * density * vel * vel) * (lk * beam) * (cf);
      float dragMag = pressureDrag + viscousDrag;

          int xtext=500;
    int ytext=216;
    int xspace=200;
    
//        text("dragmag=" + dragMag , i*xspace+xtext,ytext); ytext+=16;


      // Drag Froce
      dragMag = pressureDrag + viscousDrag;
//      dragMag = (0.5 * density * vel * vel) * (edgeLength * beam) * dragDot * 0.45;
      Vec2 dragDir = new Vec2(velDir).negate();
      Vec2 dragForce = new Vec2(dragDir).mulLocal(dragMag);
      body.applyForce( dragForce, midPoint );
      
      
      
      // Lift Force
//      float liftDot = Vec2.dot(edge, velDir); // edge^ dot velDir^ p
      float liftMag =  (0.5 * density * vel * vel) * (edgeLength * beam) * (dragDot * liftDot * 1.4144) * 1.17;
      Vec2 liftDir = Vec2.cross(1, velDir); //gets perpendicular vector (inboard pointing normal - assuming ccw body winding)
      Vec2 liftForce = new Vec2(liftDir).mulLocal(liftMag);
      body.applyForce( liftForce, midPoint );
      
      
      
      // Lifting Force in the style of Savitsky
//      float b = 2.0; // (beam - meters)
//      float beta = 0.0 * (PI/180.0); // deadrise - degrees
      
//      float cv = vel/sqrt((gravity.length()*beam)); 
//      float tau = acos(liftDot);
//      float lk = edgeLength;
//      float lambda = (lk/beam)-(1.0/(2*PI))*(tan(beta)/tan(tau+Settings.EPSILON));
      
/*
      if(lambda <= 0 || lambda > lk)
        continue;
*/    
//      float cl0 = pow(tau, 1.1)*(0.0120 * pow(lambda, 0.5) + 0.0055 * pow(lambda, 5.0/2.0) / pow(cv, 2.0));
//      float clBeta = cl0 - 0.0065 * (beta * (180.0/PI)) * pow(cl0, 0.60);
      
      // liftMag = 0.5*density*vel*vel*cl0*beam*beam;
      // liftForce = new Vec2(liftDir).mulLocal(liftMag);
      
//    int xtext=500;
//    int ytext=216;
//    int xspace=200;
    
//        text("Index=" + i , i*xspace+xtext,ytext); ytext+=16;
//        text("Total Verts" + submergedPortion.result.size() , i*xspace+xtext,ytext); ytext+=16;
//                text("beam=" + beam , i*xspace+xtext,ytext); ytext+=16;
//    text("Cv=" + cv , i*xspace+xtext,ytext); ytext+=16;
 //   text("tau=" + tau * (180.0/PI) , i*xspace+xtext,ytext); ytext+=16;
//    text("edgeLength=" + edgeLength,i*xspace+xtext,ytext); ytext+=16;
//    text("liftDot=" + liftDot, i*xspace+xtext,ytext); ytext+=16;




//    text("lk=" + lk , i*xspace+xtext,ytext); ytext+=16;
//    text("lambda=" + lambda , i*xspace+xtext,ytext); ytext+=16;
//    text("cl0=" + cl0 , i*xspace+xtext,ytext); ytext+=16;
//    text("clBeta=" + clBeta , i*xspace+xtext,ytext); ytext+=16;
//    text("liftMag=" + liftMag , i*xspace+xtext,ytext); ytext+=16;
//        text("density=" + density , i*xspace+xtext,ytext); ytext+=16;
//                text("dragdot=" + dragDot , i*xspace+xtext,ytext); ytext+=16;
//                text("maxPolygon=" + Settings.maxPolygonVertices , i*xspace+xtext,ytext); ytext+=16;


//      body.applyForce( liftForce, midPoint );
      

      
    } 



//  Statistical stuff to screen - can be commented out
fill(0);
/*
  text("Local Boat CG:" + cg,300,64);
  text("World Boat CG:" + cgWorld,300,80);
  text("World Body CG:" + body.getWorldCenter(),300,96);
  text("Trimming Torque:" + arm * buoyancyForce.length(),300,112);
  text("Trimming Angle (deg):" + body.getAngle()*(180/PI),300,128);
  text("Speed (m/s):" + new Vec2(body.getLinearVelocity()).length(), 300,144);
*/



  text("Water Area=" + submergedPortion.resultCenterArea.area,12,128+32);
  text("Water Center=" + box2d.coordWorldToPixels(submergedPortion.resultCenterArea.centroid),12,128+64);
  for(int i=0; i < submergedPortion.resultAsArray.length; ++i)
    text("    " + box2d.coordWorldToPixels(submergedPortion.resultAsArray[i]) ,12,128+64+32+i*16);

//  Draw submerged portion outline
  stroke(0,255,0);
  strokeWeight(1);
  fill(0);

  beginShape();    
  for(int i=0; i < submergedPortion.result.size(); ++i){
    Vec2 v = box2d.coordWorldToPixels(submergedPortion.resultAsArray[i]);
    vertex(v.x, v.y);
  }
  endShape(CLOSE);
}

/********************** END BUOYANCY *************************************/



  // Draw the shape associated with each fixture attached to the body - as they are in the world.
  void display() {
    fill(col,128);
    stroke(0);
    strokeWeight(1);

    Fixture fl = body.getFixtureList();
    
    while(fl != null){
      int shapeVertexCount = ((PolygonShape)fl.getShape()).getVertexCount();
      Vec2 shapeVertices[] = ((PolygonShape)fl.getShape()).getVertices();
      
      beginShape();
      for(int i=0; i< shapeVertexCount; ++i) {
        Vec2 v = box2d.coordWorldToPixels(fl.getBody().getWorldPoint(shapeVertices[i]));
        vertex(v.x, v.y);
        }
      endShape(CLOSE);
    
      fl=fl.m_next;
    }
  }
  
  
}