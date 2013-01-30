require 'nmatrix'

#extend NMATRIX::YaleFunctions

require 'pry'

@rand = NMatrix.random(4)
@yale = NMatrix.new(:yale, [4,9], 2, :int32)
0.upto(@yale.rows-1) do |row|
  0.upto(@yale.cols-1) do |column|
    @yale[row,column] = rand 100
  end
end

@list = NMatrix.new(:list, [2,2,3,4,5],  :int32)

@rand.pp
@yale.pp
puts "LIST: #{@list.pp}"

@yale * 2.0

#TODO add a #fill method to the NMatrix (for filling a testing array)

