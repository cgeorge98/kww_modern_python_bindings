#!/usr/bin/env ruby

# Auxiliary script for the parametrization of libkww:
# - read a table containing fluctuating limits omega(beta),
# - reverse the table so that beta is decreasing,
# - keep only omega that decrease,
# - so that the resulting table omega(beta) is monotonic.
# Next, we fit a simple curve to the resulting omega(beta),
# and use it as hard-coded functions kwwc_lim_hig, kwws_lim_hig.

unless ARGV.size==1
    puts "usage:\n"
    puts "   kww_lowestlims <file>"
    puts "input file:"
    puts "   must contain columns beta, w_L, w_H"
    puts "output written to stdout:"
    puts "   columns beta, w_H_enveloppe"
    exit
end

a= []
f=File.open( ARGV[0] )
f.each_line do |lin|
    a.push( lin.split.collect{ |w| w.to_f } )
end

a.reverse!

ymax = nil
a.each do |x,dummy,y|
    if !ymax || y<ymax
        ymax = y
        puts "#{x} #{y}"
    end
end
