require 'rspec'
$:.unshift( File.dirname(__FILE__))
require 'savgol'

describe SVDMatrix do 
  it "psuedoinverse of the psuedoinverse is the original matrix" do 
    a = SVDMatrix.new([2,7], [-1, 0, 1, 2, 3, 4, -7, -2, 0, 1, 1, 1, 2], :float32)
    a.pinv.pinv.should eq(a)
  end



end


describe Savgol do

  describe 'padding the ends' do

    subject do
      ar = [1, 2, 3, 4, -7, -2, 0, 1, 1]
      ar.extend(Savgol)
      ar
    end

    it 'pads with the reverse geometrically inverted sequence' do
      subject.sg_pad_ends(2).should == [-1, 0, 1, 2, 3, 4, -7, -2, 0, 1, 1, 1, 2]
      subject.sg_pad_ends(3).should == [-2, -1, 0, 1, 2, 3, 4, -7, -2, 0, 1, 1, 1, 2, 4]
    end
  end

  describe 'smoothing a signal' do
    subject do
      ar = [1, 2, 3, 4, 3.5, 5, 3, 2.2, 3, 0, -1, 2, 0, -2, -5, -8, -7, -2, 0, 1, 1]
      ar.extend(Savgol)
      ar
    end

    it "works for the simple case" do
      numpy_savgol_output = [1.0, 2.0, 3.12857143, 3.57142857, 4.27142857, 4.12571429, 3.36857143, 2.69714286, 2.04, 0.32571429, -0.05714286, 0.8, 0.51428571, -2.17142857, -5.25714286, -7.65714286, -6.4, -2.77142857, 0.17142857, 0.91428571, 1.0]
      sg = subject.savgol(5,3)
      sg.size.should == numpy_savgol_output.size
      sg.zip(numpy_savgol_output) do |sgv, numpy_sgv|
        sgv.should be_within(0.000001).of(numpy_sgv)
      end
    end
  end
end



