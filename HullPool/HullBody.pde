class HullBody {
  
  ArrayList<BuoyantBody> buoyantBodies;
  
  // Constructor
  HullBody(ArrayList<BuoyantBody> bb){
    // make a local copy of the source BuoyantBodies
    buoyantBodies = new ArrayList<BuoyantBody>(bb);
/*
    // Copy each buoyant body 
    for(int i=0; i < bb.size(); ++i) {
      buoyantBodies.add(bb.get(i));
    }
*/
    //assert(buoyantBodies.size() >= 2);

    if(buoyantBodies.size() >= 2)
    // Join each buoyant body with a welded joint    
    for(int i=1; i < buoyantBodies.size(); ++i) {
      
      WeldJointDef wjd = new WeldJointDef();
      wjd.bodyA = buoyantBodies.get(0).body;  // All bodies are joined to the 'root' 
      wjd.bodyB = buoyantBodies.get(i-0).body;  // body to minimize joint springiness stackup.
      wjd.localAnchorA.set( new Vec2(0,0) /* wjd.bodyA.getLocalCenter() */);
      wjd.localAnchorB.set(new Vec2(0,0) /* wjd.bodyB.getLocalCenter() */);
      wjd.referenceAngle = 0.0;

      wjd.collideConnected=false;
      // These properties affect how springy the joint is 
      wjd.frequencyHz = 0;  // Try a value less than 5 (0 for no elasticity)
      wjd.dampingRatio = 1; // Ranges between 0 and 1 (1 for no springiness)

      // Make the joint
      WeldJoint wj = (WeldJoint) box2d.world.createJoint(wjd);

      
//      RevoluteJointDef wjd = new RevoluteJointDef();
//      wjd.collideConnected=false;
//      wjd.bodyA = buoyantBodies.get(0).body;  // All bodies are joined to the 'root' 
//      wjd.bodyB = buoyantBodies.get(i-0).body;  // body to minimize joint springiness stackup.
//      wjd.localAnchorA.set( new Vec2(0,0) /* wjd.bodyA.getLocalCenter() */);
//      wjd.localAnchorB.set(new Vec2(0,0) /* wjd.bodyB.getLocalCenter() */);
//      wjd.referenceAngle = 0.0;
      
//      wjd.enableLimit=true;
//      wjd.lowerAngle=0;
//      wjd.upperAngle=0.0000001;

//      PrismaticJointDef wjd = new PrismaticJointDef();
//      wjd.collideConnected=false;
//      wjd.bodyA = buoyantBodies.get(0).body;  // All bodies are joined to the 'root' 
//      wjd.bodyB = buoyantBodies.get(i-0).body;  // body to minimize joint springiness stackup.
//      wjd.localAnchorA.set( new Vec2(0,0) /* wjd.bodyA.getLocalCenter() */);
//      wjd.localAnchorB.set(new Vec2(0,0) /* wjd.bodyB.getLocalCenter() */);

//      wjd.referenceAngle = 0.0;
      
//      wjd.enableLimit=true;
//      wjd.lowerTranslation=0.0;
//      wjd.upperTranslation=0.0;

      

      // Make the joint
      // WeldJoint wj = (WeldJoint) box2d.world.createJoint(wjd);

      
    }

  }
  
  
  
    // Draw the shape associated with each BuoyantBody.
  void processAndDisplay() {

    // Must display first
    for(int i=0; i < buoyantBodies.size(); ++i) {
      buoyantBodies.get(i).display();
    }
    
    // Then process buoyancy
    for(int i=0; i < buoyantBodies.size(); ++i) {
      buoyantBodies.get(i).processBuoyancy();
    }

    
  }
  
  
}
