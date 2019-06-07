
// The code in this tab is used to create all the other
// controls needed to configure the slider control.


GSlider sldrSigWaveHt, sldrTargetVelocity;
GLabel lblSigWaveHt, lblTargetVelocity;

public void makeSliderConfigControls() {
  // Create colour scheme selectors
  int x,y;
  
  x = width - 180; 
  y = 100;

  // Significant Wave Ht Slider & Lavel
  sldrSigWaveHt = new GSlider(this, x, y + (34*1), 80, 40, 12);
  sldrSigWaveHt.setLimits(0.0f, 0.0f, 6.0f); // Default Val, Lo Val, Hi Val
  sldrSigWaveHt.setShowValue(true);
 
  lblSigWaveHt = new GLabel(this, x + 82, y + (34*1), 80, 40, "(h) Sig. Wave Ht - m");
  lblSigWaveHt.setTextAlign(GAlign.LEFT, null);
  lblSigWaveHt.setTextItalic();

  // Target Velocity Slider  & Label
  sldrTargetVelocity = new GSlider(this, x, y + (34*2), 80, 40, 12);
  sldrTargetVelocity.setLimits(0.001f, 0.001f, 30.0f); // Default Val, Lo Val, Hi Val
  sldrTargetVelocity.setShowValue(true);

  lblTargetVelocity = new GLabel(this, x + 82, y + (34*2), 80, 40, "(V) Target Velocity - m/s");
  lblTargetVelocity.setTextAlign(GAlign.LEFT, null);
  lblTargetVelocity.setTextItalic();

}