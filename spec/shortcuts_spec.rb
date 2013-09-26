# = NMatrix
#
# A linear algebra library for scientific computation in Ruby.
# NMatrix is part of SciRuby.
#
# NMatrix was originally inspired by and derived from NArray, by
# Masahiro Tanaka: http://narray.rubyforge.org
#
# == Copyright Information
#
# SciRuby is Copyright (c) 2010 - 2012, Ruby Science Foundation
# NMatrix is Copyright (c) 2012, Ruby Science Foundation
#
# Please see LICENSE.txt for additional copyright notices.
#
# == Contributing
#
# By contributing source code to SciRuby, you agree to be bound by
# our Contributor Agreement:
#
# * https://github.com/SciRuby/sciruby/wiki/Contributor-Agreement
#
# == shortcuts_spec.rb
#
# Specs for the shortcuts used in NMatrix and in NVector.
#

# Can we use require_relative here instead?
require File.join(File.dirname(__FILE__), "spec_helper.rb")
require 'pry'

describe NMatrix do

  it "zeros() creates a matrix of zeros" do
    m = NMatrix.zeros(3)
    n = NMatrix.new([3, 3], 0)

    m.should.eql? n
  end

  it "ones() creates a matrix of ones" do
    m = NMatrix.ones(3)
    n = NMatrix.new([3, 3], 1)

    m.should.eql? n
  end

  it "eye() creates an identity matrix" do
    m = NMatrix.eye(3)
    identity3 = NMatrix.new([3, 3], [1, 0, 0, 0, 1, 0, 0, 0, 1])

    m.should.eql? identity3
  end

  it "diag() creates a matrix with pre-supplied diagonal" do
    arr = [1,2,3,4]
    m = NMatrix.diag(arr)
    m.is_a?(NMatrix).should be_true
  end

  it "diagonals() contains the seeded values on the diagonal" do
    arr = [1,2,3,4]
    m = NMatrix.diagonals(arr)
    m[0,0].should eq(arr[0])
    m[1,1].should eq(arr[1])
    m[2,2].should eq(arr[2])
    m[3,3].should eq(arr[3])
  end
  it "diagonals() allows you to set the dtype per the documentation" do 
    arr = [1,2,3,4]
    m = NMatrix.diagonals(arr, dtype: :float32)
    m.dtype.should eq(:float32)
    m[0,0].should eq(arr[0])
    m[1,1].should eq(arr[1])
    m[2,2].should eq(arr[2])
    m[3,3].should eq(arr[3])
  end

  it "random() creates a matrix of random numbers" do
    m = NMatrix.random(2)

    m.stype.should == :dense
    m.dtype.should == :float64
  end

  it "random() only accepts an integer or an array as dimension" do
    m = NMatrix.random([2, 2])

    m.stype.should == :dense
    m.dtype.should == :float64

    expect { NMatrix.random(2.0) }.to raise_error
    expect { NMatrix.random("not an array or integer") }.to raise_error
  end

  it "seq() creates a matrix of integers, sequentially" do
    m = NMatrix.seq(2) # 2x2 matrix.
    value = 0

    2.times do |i|
      2.times do |j|
        m[i,j].should == value
        value += 1
      end
    end
  end


  it "indgen() creates a matrix of integers as well as seq()" do
    m = NMatrix.indgen(2) # 2x2 matrix.
    value = 0

    2.times do |i|
      2.times do |j|
        m[i, j].should == value
        value += 1
      end
    end
  end

  it "findgen creates a matrix of floats, sequentially" do
    m = NMatrix.findgen(2) # 2x2 matrix.
    value = 0

    2.times do |i|
      2.times do |j|
        (m[i, j]/10).should be_within(Float::EPSILON).of(value.to_f/10)
        value += 1
      end
    end
  end

  it "bindgen() creates a matrix of bytes" do
    m = NMatrix.bindgen(2) # 2x2 matrix.
    value = 0

    2.times do |i|
      2.times do |j|
        m[i, j].should == value
        value += 1
      end
    end
  end

  it "cindgen() creates a matrix of complexes" do
    m = NMatrix.cindgen(2) # 2x2 matrix.
    value = 0

    2.times do |i|
      2.times do |j|
        m[i, j].real.should be_within(Float::EPSILON).of(value)
        m[i, j].imag.should be_within(Float::EPSILON).of(0.0)
        value += 1
      end
    end
  end

  it "column() returns a NMatrix" do
    m = NMatrix.random(3)

    m.column(2).is_a?(NMatrix).should be_true
  end

  it "row() returns a NMatrix" do
    m = NMatrix.random(3)

    m.row(2).is_a?(NMatrix).should be_true
  end

  it "diagonals() creates an NMatrix" do
    arr = [1,2,3,4]
    m = NMatrix.diagonals(arr)
    m.is_a?(NMatrix).should be_true
  end

  it "diagonals() contains the seeded values on the diagonal" do
    arr = [1,2,3,4]
    m = NMatrix.diagonals(arr)
    m[0,0].should eq(arr[0])
    m[1,1].should eq(arr[1])
    m[2,2].should eq(arr[2])
    m[3,3].should eq(arr[3])
  end

  context "_like constructors" do
    before :each do
      STDERR.puts "starting GC"
      GC.start
      STDERR.puts "GC finished"
      @nm_1d = NMatrix[5.0,0.0,1.0,2.0,3.0]
      STDERR.puts "@nm_1d"
      @nm_2d = NMatrix[[0.0,1.0],[2.0,3.0]]
      STDERR.puts "@nm_2d"
    end

    it "should create an nmatrix of ones with dimensions and type the same as its argument" do
      NMatrix.ones_like(@nm_1d).should eq NMatrix[1.0, 1.0, 1.0, 1.0, 1.0]
      NMatrix.ones_like(@nm_2d).should eq NMatrix[[1.0, 1.0], [1.0, 1.0]]
    end

    it "should create an nmatrix of zeros with dimensions and type the same as its argument" do
      STDERR.puts "A"
      NMatrix.zeros_like(@nm_1d).should eq NMatrix[0.0, 0.0, 0.0, 0.0, 0.0]
      STDERR.puts "B"
      NMatrix.zeros_like(@nm_2d).should eq NMatrix[[0.0, 0.0], [0.0, 0.0]]
      STDERR.puts "C"
    end
  end

end

describe "NVector" do

  it "zeros() creates a vector of zeros" do
    v = NVector.zeros(4)

    4.times do |i|
      v[i].should == 0
    end
  end

  it "ones() creates a vector of ones" do
    v = NVector.ones(3)

    3.times do |i|
      v[i].should == 1
    end
  end

  it "random() creates a vector of random numbers" do
    v = NVector.random(4)
    v.dtype.should == :float64
    v.stype.should == :dense
  end

  it "seq() creates a vector of integers, sequentially" do
    v = NVector.seq(7)
    v.should == NMatrix.new([7,1], [0, 1, 2, 3, 4, 5, 6])
  end

  it "seq() only accepts integers as dimension" do
    expect { NVector.seq(3) }.to_not raise_error

    expect { NVector.seq([1, 3]) }.to raise_error
    expect { NVector.seq(:wtf) }.to raise_error
  end

  it "indgen() creates a vector of integers as well as seq()" do
    v = NVector.indgen(7)
    v.should == NMatrix.new([7,1], [0, 1, 2, 3, 4, 5, 6])
  end

  it "findgen creates a vector of floats, sequentially" do
    v = NVector.findgen(2)
    v.should == NMatrix.new([2,1], [0.0, 1.0])
  end

  it "bindgen() creates a vector of bytes, sequentially" do
    v = NVector.bindgen(4)
    v.should == NMatrix.new([4,1], [0, 1, 2, 3], dtype: :byte)
  end

  it "cindgen() creates a vector of complexes, sequentially" do
    v = NVector.cindgen(2)
    v.should == NMatrix.new([2,1], [Complex(0.0, 0.0), Complex(1.0, 0.0)], dtype: :complex64)
  end

  it "linspace() creates a vector with n values equally spaced between a and b" do
    v = NVector.linspace(0, 2, 5)
    v.should == NMatrix.new([5,1], [0.0, 0.5, 1.0, 1.5, 2.0])
  end

  it "logspace() creates a vector with n values logarithmically spaced between decades 10^a and 10^b" do
    v = NVector.logspace(0, 3, 4)
    v.should == NMatrix.new([4,1], [1.0, 10.0, 100.0, 1000.0])
  end
end

describe "Inline constructor" do

  it "creates a NMatrix with the given values" do
    m = NMatrix.new([2, 2], [1, 4, 6, 7])
    n = NMatrix[[1, 4], [6, 7]]

    m.should.eql? n
  end
end
