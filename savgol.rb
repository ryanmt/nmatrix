$:.unshift "lib"
require 'nmatrix'

require 'pry'

JTP_Sigma = true
JTP_SVD = true

class SVDMatrix < NMatrix
  def pinv
    (s_temp,u,v) = NMatrix::LAPACK.svd(self, :all) 
    # This differs from the example, as we run into size differences with the matrix outputs from SVD,
    # as U is supposed to be MxM and V* is NxN.  S(sigma) is returned from LAPACK as just diagonal elements, 
    # so I'm filling it into the MxN shape for the later processing
    s = NMatrix.new( [u.shape.first, v.shape.first], 0, u.dtype)
    p s_temp.to_a.map(&:first)
    if JTP_Sigma or JTP_SVD
      s_temp = [1.3733482360839844, 6.092119216918945, 11.788289070129395, 1.017959475517273 ]
      s_temp.each_with_index do |n,i|
        s[i,i] = n
      end
      if JTP_SVD
        u = NMatrix.new(u.shape, [0.17072263360023499, 0.6756982803344727, 0.6994683742523193, -0.10365278273820877, -0.4787352979183197, 0.20552638173103333, 0.10365264117717743, 0.6994684338569641, -0.6952213048934937, 0.048802461475133896, 4.372066442215328e-08, 1.4553515725879151e-08, -0.4787352681159973, 0.2055264115333557, -0.10365280508995056, -0.6994684338569641, 0.17072266340255737, 0.6756982207298279, -0.6994683146476746, 0.10365286469459534] , u.dtype)
        v = NMatrix.new(v.shape, [-0.9547808170318604, 0.29731065034866333, 0.0, -0.0, 0.0, 0.0, -0.2549290359020233, -0.9669595956802368, 0.29731065034866333, 0.9547808170318604, -3.952681382202172e-08, 4.46132837339519e-08, 0.0, 0.0, -0.9669598340988159, 0.25492921471595764] , v.dtype)
      end
    else 
      s_temp.col(0).to_a.each_with_index do |n, i| 
        s[i,i] = n.first 
      end
    end
    puts "S: "
    pp s
    puts "A: #{self.shape}, S: #{s.shape}, U: #{u.shape}, V: #{v.shape}"
    puts "U: "
    pp u
    puts "V: "
    pp v

    # PINV... if A = U . S . (V*) then A+ = V . S+ . (U*).  
    # Take the reciprocal of each non-zero diagonal element of S (which is all we have, so just the reciprocal),
    # then do the reversal and processing of the other elements
    # SOURCE: http://en.wikipedia.org/wiki/Moore%E2%80%93Penrose_pseudoinverse#Singular_value_decomposition_.28SVD.29
    if u.dtype == :complex64 or u.dtype == :complex128
      v.conjugate_transpose.dot(s.inverse.tranpose).dot(w.conjugate_transpose)
    else
      # "Finding the conjugate transpose of a matrix A with real entries reduces to finding the transpose of A, as the conjugate of a real number is the number itself." http://en.wikipedia.org/wiki/Conjugate_transpose
      # v.transpose.dot(s.transpose.dot(s).inverse).dot(u.transpose) # Fails for some reason...
      #v.transpose.dot(s.inverse.transpose).dot(u.transpose)
     NMatrix::BLAS.gemm(v.transpose, NMatrix::BLAS.gemm(s.inverse.transpose, u.transpose))
    end
  end
end


module Savgol
  def sg_check_args(window_size, order)
    # checks that the window size is a positive and odd integer
    if !window_size.is_a?(Integer) || window_size < 1 || window_size % 2 != 1
      raise ArgumentError, "window_size must be a positive odd integer"
    end
    if !order.is_a?(Integer) || order < 0
      raise ArgumentError, "order must be an integer >= 0"
    end
    if window_size < order + 2
      raise ArgumentError, "window_size is too small for the polynomial order"
    end
  end

  def savgol(window_size, order, deriv = 0, check_args=true)
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
    start.map! {|v| self[0] - (v-self[0]).abs }

    fin = self[(-half_window-1)...-1]
    fin.reverse!
    fin.map! {|v| self[-1] + (v-self[-1]).abs }
    start.push(*self, *fin)
  end

  # returns an object that will convolve with the padded array
  def sg_weights(half_window, order, deriv = 0)
    mat = SVDMatrix.new([order + 2, half_window*2], (-half_window..half_window).map {|k| (0..order).map {|i| k**i }}.flatten, :float32)
    #mat.pinv.row(deriv).to_a
    resp = mat.pinv
    pp resp
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
