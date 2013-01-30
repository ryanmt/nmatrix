require 'spec_helper.rb'
#  Types
DENSE = NMatrix.random(4)
YALE = NMatrix.new(:yale, [4,9], 2, :int32)
0.upto(YALE.rows-1) do |row|
  0.upto(YALE.cols-1) do |column|
    YALE[row,column] = rand 100
  end
end
LIST = NMatrix.new(:list, [2,2,3,4,5],  :int32)

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
  test_nms = {dense: DENSE, yale: YALE, list: LIST }
  [:dense, :yale, :list].each do |dtype|
    describe ":#{dtype.upcase}" do 
      it "takes a scalar" do 
        @test = test_nms[dtype]
        @prod = @test * 2.0
        @prod[1,1].should == (@test[1,1] * 2.0)
      end
      it "takes a matrix (C code)" do 
      end
    end
  end
end
