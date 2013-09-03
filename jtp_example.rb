# do what the heck you want to public license (see doc end)
 
# gem install ruby-svg   # provides SVDMatrix
require 'ruby-svd'

NMatrix_SVD = false

require 'pry'
class Matrix
  def pretty_print
    i = 0
    j = 0
    print "[\n  ["
    self.each do |number|
      print "  " + "%.9f" % number.to_s
      i+= 1
      if i == self.column_size
        j+= 1
        if j == self.row_size
          print "]\n]\n"
        else
          print "]\n"
          print "  [" 
        end
        i = 0
      else
        print ", "
      end
    end
  end
  def shape
    [self.row_size, self.column_size]
  end
end 
class SVDMatrix < Matrix
  def self.[](*rows)
    mat = self.new(rows.size,rows.first.size)
    rows.each_with_index {|row,i| mat.set_row(i, row) }
    mat
  end
 
  # Moore-Penrose psuedo-inverse
  # SVD: A=UWV'
  # A+=V(W'W)^(−1)W'U'
  def pinv
    (u, w, v) = self.decompose
    if NMatrix_SVD
      w = Matrix[ [10.287137031555176, 0.0, 0.0, 0.0],
        [ 0.0, 8.434704780578613, 0.0,  0.0],
        [ 0.0, 0.0, 1.220271348953247,  0.0],
        [ 0.0, 0.0, 0.0, 0.7358231544494629],
        [ 0.0, 0.0, 0.0,  0.0]]
        u = Matrix[[-0.012036953121423721,   0.215118408203125, -0.13782189786434174,   0.8375576138496399,  0.48278582096099854],
          [ -0.18088601529598236,  0.0964784175157547,  -0.5232441425323486,   0.3249719440937042,  -0.7606456875801086],
          [  -0.5396060347557068,  0.5359153747558594,  -0.4121575951576233, -0.39301204681396484,   0.3119094967842102],
          [  -0.6421598792076111, 0.15470080077648163,   0.6912844181060791,   0.1922229677438736, -0.22107598185539246],
  [  -0.5134116411209106, -0.7957878708839417, -0.24387042224407196, 0.038505882024765015,   0.2053646445274353]]
        v = Matrix[
          [-0.7009949684143066,  -0.6908645629882812, 0.16057008504867554, 0.07436076551675797],
          [ 0.1168968677520752,  0.13344630599021912,  0.8970689177513123,  0.4047156274318695],
          [0.11495141685009003, -0.11368716508150101, -0.4042222201824188,  0.9002586603164673],
          [ 0.6940658092498779,  -0.7014082670211792, 0.07803323119878769, -0.1421615034341812]
        ]
    end
    puts w
    puts "S: "
    w.pretty_print
    puts "A: #{self.shape}, S: #{w.shape}, U: #{u.shape}, V: #{v.shape}"
    puts "U: "
    u.pretty_print
    puts "V: "
    v.pretty_print
    binding.pry
    v * (w.t*w).inverse * w.t * u.t
  end
end
 
module Savgol
  def sg_check_args(window_size, order)
    # checks that the window size is a positive odd Integer
    if !window_size.is_a?(Integer) || window_size < 1 || window_size % 2 != 1
      raise ArgumentError, "window_size size must be a positive odd integer" 
    end
    if !order.is_a?(Integer) || order < 0
      raise ArgumentError, "order must be an integer >= 0"
    end
    if window_size < order + 2
      raise ArgumentError, "window_size is too small for the polynomial order"
    end
  end
 
  def savgol(window_size, order, deriv=0, check_args=true)
    sg_check_args(window_size, order) if check_args
    half_window = (window_size -1) / 2
    weights = sg_weights(half_window, order, deriv)
    ar = sg_pad_ends(half_window)
    sg_convolve(ar, weights)
  end
 
  def sg_convolve(data, weights, mode=:valid)
    data.each_cons(weights.size).map do |ar|
      ar.zip(weights).map {|pair| pair[0] * pair[1] }.reduce(:+)
    end
  end
 
  # pads the ends with the reverse, geometric inverse sequence
  def sg_pad_ends(half_window)
    start = self[1..half_window]
    start.reverse!
    start.map! {|v| self[0] - (v - self[0]).abs }
 
    fin = self[(-half_window-1)...-1]
    fin.reverse!
    fin.map! {|v| self[-1] + (v - self[-1]).abs }
    start.push(*self, *fin)
  end
 
  # returns an object that will convolve with the padded array
  def sg_weights(half_window, order, deriv=0)
    mat = SVDMatrix[ *(-half_window..half_window).map {|k| (0..order).map {|i| k**i }} ]
    #mat.pinv.row(deriv).to_a
    resp = mat.pinv
    resp.pretty_print
    resp.row(deriv).to_a
  end
end
 
class Array
  include Savgol
end
 
if __FILE__ == $0
  ar = [1, 2, 3, 4, 3.5, 5, 3, 2.2, 3, 0, -1, 2, 0, -2, -5, -8, -7, -2, 0, 1, 1] 
  smoothed = ar.savgol(5,3)
  # smoothed => [0.9999999725324346, 2.0000000115857715, 3.1285714305655055, 3.571428673088599, 4.271428614306785, 4.1257144088160675, 3.3685715767719575, 2.6971429967311664, 2.0400000845020303, 0.3257144724996456, -0.057142738662666615, 0.7999999264680943, 0.5142855346070793, -2.1714285517550818, -5.257142779184832, -7.657142699240195, -6.40000007190543, -2.771428896122372, 0.17142820252085061, 0.9142855954984234, 1.0000000154252753]
  puts smoothed
end
