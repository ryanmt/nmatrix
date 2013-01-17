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

@rand.pp
@yale.pp


#TODO add a #fill method to the NMatrix (for filling a testing array)

# SEG FAULT?
# @yale[3,4] = rand(3)
# /home/ryanmt/.rvm/gems/ruby-1.9.2-p320/gems/pry-0.9.10/lib/pry/indent.rb:190: [BUG] Segmentation fault

