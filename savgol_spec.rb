require 'rspec'
require 'savgol'

describe SVDMatrix do 
  it "psuedoinverse of the psuedoinverse is the original matrix" do 
    
    a.pinv.pinv.should eq(a)
  end
end


