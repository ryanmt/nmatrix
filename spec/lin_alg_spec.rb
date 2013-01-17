require 'spec_helper.rb'

describe "Vector Norms" do 
  pending "p-norms"
  pending 'weighted p-norms'
end

describe "Matrix Norms" do 
  pending 'induced matrix norm'
  pending "p-norms"
  pending "weighted p-norms"
end

describe "Matrix-Scalar multiplication" do 
  it "takes a scalar" do 
    @rand = NMatrix.random(4)
    @prod = @rand * 2.0
    @prod[1,1].should == (@rand[1,1] * 2.0)
  end
  it "takes a matrix (C code)" do 
    @mat1 = NMatrix.identity(3)
    @mat2 = NMatrix.new([3,3], [1,2,3])
    @mat_result = @mat1 * @mat2
    @mat_result.pp
  end
end
