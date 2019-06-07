// An deep water ocean wave body in JBOX2D

class WaveBody {
  Body body; // keep track of the body.
  
  // We'll keep track of all of the surface points
  ArrayList<Vec2> surface;
  
  // Wave characteristic values;
  float worldGravity;

  float wavePeriod;
  float h13;
  
  float waveLength;
  float waveNumber;
  float waveSpeed;
  
  boolean waveDirection;
  
  float fetchLength;
  float fetchStep;
  
  float depth;
  
  Vec2 position;
  
  
  WaveBody() {
    position = new Vec2();
    body = null;
    surface = new ArrayList<Vec2>();

    // Default Values
    
      worldGravity = box2d.world.getGravity().length();

    
    wavePeriod=2.0; // (second)
    h13=1.0; // (meter)
    
    depth=10; // (meter)
      
    waveLength=worldGravity*wavePeriod*wavePeriod/TWO_PI; // (meter)
    waveNumber=TWO_PI/waveLength; // (rad / meter)
    waveSpeed=waveLength*(1.0/wavePeriod); // (meter / second)
    
    waveDirection=false;
      
    fetchLength=50; // (meter)
    fetchStep=0.25; // (meter)
    
    float t=0.0;
    
    AdvanceWave(t);
  }


WaveBody(float posX, float posY) {
    position = new Vec2(posX, posY);
    body = null;
    surface = new ArrayList<Vec2>();

    // Default Values
    
     worldGravity = box2d.world.getGravity().length();
    
    wavePeriod=2.0; // (second)
    h13=1.0; // (meter)
    depth=10; // (meter)
    
    waveLength=worldGravity*wavePeriod*wavePeriod/TWO_PI; // (meter)
    waveNumber=TWO_PI/waveLength; // (rad / meter)
    waveSpeed=waveLength*(1.0/wavePeriod); // (meter / second)

    waveDirection=false;
      
    fetchLength=50; // (meter)
    fetchStep=0.25; // (meter)
    
    float t=0.0;
    
    AdvanceWave(t);

}

void SetWaveParams(float _period, float _h13) {
    wavePeriod = _period; // (second)
    h13 = _h13; // (meter)
      
    waveLength=worldGravity*wavePeriod*wavePeriod/TWO_PI; // (meter)
    waveNumber=TWO_PI/waveLength; // (rad / meter)
    waveSpeed=waveLength*(1.0/wavePeriod); // (meter / second)

}




  // This function removes the particle from the box2d world
  void killBody() {
    box2d.destroyBody(body);
  }


// Advance the wavebody to time t
void AdvanceWave(float t){

    if(body != null){
      killBody();
      surface.clear();
    }

    
    // This has to go backwards (counter-clockwise winding) so that the objects  bounce off the top of the surface
    // This "edgechain" will only work in one direction!
    // 
    // Wave Surface
    for (float x = fetchLength; x >= 0; x -= fetchStep) {

      // Calculate the surface at each point x at time t
      float y=(h13/2)*cos((TWO_PI/waveLength)*(x-((waveDirection)? waveSpeed: -waveSpeed)*t));
  
      // Store the vertex in screen coordinates
      surface.add(new Vec2(x,y));
    }
    
    // Depth Box
    surface.add(new Vec2(surface.get(surface.size()-1).x, -depth));
    surface.add(new Vec2(surface.get(0).x, -depth));
//    surface.add(new Vec2(surface.get(0).x, surface.get(0).y));
    

    // Build an array of vertices in Box2D coordinates
    // from the ArrayList we made
    Vec2[] vertices = new Vec2[surface.size()];
    for (int i = 0; i < vertices.length; i++) {
      Vec2 edge = surface.get(i);
      vertices[i] = edge;
    }

    // This is what box2d uses to put the surface in its world
    ChainShape chain = new ChainShape();
    
    // Create the chain!
    chain.createLoop(vertices,vertices.length);

    // The edge chain is now attached to a body via a fixture
    BodyDef bd = new BodyDef();
    bd.type = BodyType.STATIC;
    bd.position.set(position.x,position.y);
    
    body = box2d.createBody(bd);
    
    // Attached the shape to the body using a Fixture
    FixtureDef fd = new FixtureDef();
    fd.shape=chain;
    fd.setSensor(true);  /***************** MAKE FIXTURE A SENSOR ************************/
    fd.setDensity(1024.0);
    
    body.createFixture(fd);

    body.setUserData(this);  
}




  // Draw the shape associated with each fixture attached to the body - as they are in the world.
  void display() {
    stroke(0);
    strokeWeight(1);
    //noFill();
    fill(67,197,255,64);

    Fixture fl = body.getFixtureList();
    
    while(fl != null){
      int shapeVertexCount = ((ChainShape)fl.getShape()).m_count;
      Vec2 shapeVertices[] = ((ChainShape)fl.getShape()).m_vertices;
      
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