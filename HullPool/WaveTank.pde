// A container solid for WaveBody(s) in JBOX2D

class WaveTank {
  Body body; // keep track of the body.
  
  // We'll keep track of all of the surface points
  ArrayList<Vec2> surface;
  
  // Tank characteristic values;
  float tankDepth; // Interior Depth (meters)
  float tankLength; // Interior Length (meters)
  float borderThickness; // Wall Thickness (meters) 
  Vec2 position;
  

// 
void setDefaultValues(){
    position = new Vec2();
    body = null;
    surface = new ArrayList<Vec2>();

    // Default Values
    tankDepth=12;
    tankLength=50;
    borderThickness=1;
}

  
// Create a default tank at (0,0)  
  WaveTank() {
    setDefaultValues();

    CreateTank();    
  }



// Create a default tank at the given position
WaveTank(float posX, float posY) {
    setDefaultValues();

    position.set(posX, posY);

    CreateTank();    
}



// Create a fully specified tank
WaveTank(float posX, float posY, float _depth, float _length, float _borderThickness) {
    setDefaultValues();

    position.set(posX, posY);
    tankDepth=_depth;
    tankLength=_length;
    borderThickness=_borderThickness;

    CreateTank();    
}



// This function removes the body from the box2d world
  void killBody() {
    box2d.destroyBody(body);
  }


// Create Tank
void CreateTank(){

    if(body != null){
      killBody();
      surface.clear();
    }


   
    // Tank Box
    surface.add(new Vec2(0,0));
    surface.add(new Vec2(tankLength+borderThickness,0));
    surface.add(new Vec2(tankLength+borderThickness,tankDepth+borderThickness));
    surface.add(new Vec2(tankLength,tankDepth+borderThickness));
    
    surface.add(new Vec2(tankLength, borderThickness));
    surface.add(new Vec2(borderThickness,borderThickness));
    surface.add(new Vec2(borderThickness,tankDepth+borderThickness));
    surface.add(new Vec2(0,tankDepth+borderThickness));
    surface.add(new Vec2(0,0));
    

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
    chain.createChain(vertices,vertices.length);

    // The edge chain is now attached to a body via a fixture
    BodyDef bd = new BodyDef();
    bd.type = BodyType.STATIC;
    bd.position.set(position.x,position.y);
    
    body = box2d.createBody(bd);
    
    // Attached the shape to the body using a Fixture
    FixtureDef fd = new FixtureDef();
    fd.shape=chain;
    
    body.createFixture(fd);

    body.setUserData(this);  
}




 // Draw the shape associated with each fixture attached to the body - as they are in the world.
  void display() {
    stroke(0);
    strokeWeight(1);
    fill(0,32);

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
