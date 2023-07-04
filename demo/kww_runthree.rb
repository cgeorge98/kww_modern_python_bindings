#!/usr/bin/env ruby

# run all three algorithms (l|m|h)

def callkww( argstr )
    system( "runkww 0 #{argstr} > /ram/b 2> /ram/c" )
    c = `cat /ram/c`
    c != "" and
        raise "invalid result:\n#{c}"
    return `cat /ram/b`.split
end

unless ARGV.size==3
    puts "usage:\n"
    puts "   #{$0} s|c <b> <w>"
    puts "with arguments:"
    puts "   flag1: s: sin transform kwws"
    puts "          c: cos transform kwwc"
    puts "   b: stretching exponent (between 0.1 and 2)"
    puts "   w: omega"
    puts "return values:"
    puts "   value  1:   regular algorithm"
    puts "   values 2-4: the computed function for algorithms l,m,h"
    puts "   values 5-7: number of summed terms for algorithms l,m,h"
    exit
end

dir = ARGV[0]
b   = ARGV[1]
w   = ARGV[2]

ra = callkww( "#{dir} a #{b} #{w}" )
rl = callkww( "#{dir} l #{b} #{w}" )
rm = callkww( "#{dir} m #{b} #{w}" )
rh = callkww( "#{dir} h #{b} #{w}" )

puts "#{ra[1]}  #{rl[0]} #{rm[0]} #{rh[0]}   #{rl[2]} #{rm[2]} #{rh[2]}"
