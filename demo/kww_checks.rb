#!/usr/bin/env ruby1.9

include Math

##############################################################################
##
## auxiliary scripts for the library libkww
##
## Author::    Joachim Wuttke <j.wuttke@fz-juelich.de>
## Copyright:: 2012
## Licence::   GPL v3 or higher
##
## This script checks for the continuity of kwws or kwwc as function of beta
## or omega, concentrating on points where the used algorithm or the number
## of summed terms changes.
##
##############################################################################


#! Call runkww with given arguments, check for errors, return kwws or kwwc.

def kww( fou, alg, sca, p, v )
    argstr = case sca
             when "b" then "#{v} #{p}"
             when "w" then "#{p} #{v}"
             else raise "BUG: invalid sca"
             end
    cmd = "runkww 0 #{fou} #{alg} #{argstr} > /ram/b 2> /ram/c"
    system( cmd )
    c = `cat /ram/c`
    c != "" and
        raise "command '#{cmd}' yields invalid result:\n#{c}"
    ret = `cat /ram/b`.split()
    # puts "#{cmd} -> #{ret}"
    v = ret[0].to_f
    if v==0
        puts "#{v} -> #{ret}"
        exit
    end
    return [ ret[0].to_f, ret[1].to_i, ret[2].to_i ]
end


#! One scan.

#  Input:
#    eps:      floating-point step
#    fou:      kind of Fourier transform
#    sca:      the scan variable
#    slo:      the fixed parameter
#    p:        the value of the fixed parameter
#    scanlims: limits for the scan variable
#  No output, except for progress report and error messages

def scan( eps, fou, sca, slo, p, scanlims )
    puts "scan in #{sca} with #{slo}=#{p}"
    v1 = scanlims[0] # the scan variable
    y1, a1, n1 = kww( fou, :a, sca, p, v1 )
    while true
        v0, y0, a0, n0 = [ v1, y1, a1, n1 ]
        v0 >= scanlims[1] and return
        v1 = [ 2 * v0, scanlims[1] ].min
        y1, a1, n1 = kww( fou, :a, sca, p, v1 )
        if a1!=a0 || n1!=n0 # there is at least one change between v0 and v1
            # bisection to localize earliest change at v2
            i = 0
            while true
                v2 = (v0 + v1)/2
                i += 1
                if i>100
                    puts "bisection failed"
                    break
                end
                y2, a2, n2 = kww( fou, :a, sca, p, v2 )
                v1-v0 < eps*v2 and break
                if a2==a0 && n2==n0
                    v0, y0, a0, n0 = [ v2, y2, a2, n2 ]
                else
                    v1, y1, a1, n1 = [ v2, y2, a2, n2 ]
                end
            end
            # puts( "change at %13.8g: #{a0}:#{n0} -> #{a1}:#{n1}" % v2 )

            # points on both sides of v2
            y3, a3, n3 = kww( fou, :a, sca, p, v3 = v2*(1-6*eps) )
            y4, a4, n4 = kww( fou, :a, sca, p, v4 = v2*(1-2*eps) )
            y5, a5, n5 = kww( fou, :a, sca, p, v5 = v2*(1+2*eps) )
            y6, a6, n6 = kww( fou, :a, sca, p, v6 = v2*(1+6*eps) )
            a3==a4 && n3==n4 or
                puts "unexpected change at -6eps..-2eps"
            a5==a6 && n5==n6 or
                puts "unexpected change at +2eps..+6eps"
            a4==a5 && n4==n5 and
                puts "change at -2eps..2eps not reproducible"

            # check continuity of y
            s3 = y4-y3
            s4 = y5-y4
            s5 = y6-y5
            info = "%13.8g: #{a4}:#{n4} -> #{a5}:#{n5} steps %10.5g %10.5g %10.5g" % [ v2, s3, s4, s5 ]
            if ( s3<0 && s5<0 && s4>0 ||
                 s3>0 && s5>0 && s4<0 )
                puts "non-monoton at #{info}"
            elsif( s4.abs > s3.abs + s5.abs + eps*(y4+y5).abs )
                puts "big step    at #{info}"
            end

            # compare series expansion wit integration
            [ [v4,y4,a4,n4], [v5,y5,a5,n5] ].each do |v,y,a,n|
                a<10 and next
                yy, aa, nn = kww( fou, :m, sca, p, v )
                if yy<0
                    puts "integration at #{v} failed: #{yy}"
                elsif (yy-y).abs>eps*y
                    puts "discrepancy at #{v} (#{a}:#{n} vs #{aa}:#{nn}):" +
                        " y=#{y}, rel err #{(yy-y)/y}"
                end
            end
        end
    end
end


#! Main program.

unless ARGV.size==5
    puts "usage:"
    puts "   kwwtest1 <eps> c|s|p b|w r|s|x <n>|<val>"
    puts "with arguments:"
    puts "   <eps>: required floating-point resolution"
    puts "   flag1: c: cos transform"
    puts "          s: sin transform"
    puts "          p: primitive of cos transform"
    puts "   flag2: b: scan beta"
    puts "          w: scan omega"
    puts "   flag3: r: non-scan parameter at random"
    puts "          s: non-scan parameter in regular slices"
    puts "          x: non-scan parameter fixed"
    puts "   <n>:   number of scans (if flag3=r|s)"
    puts "   <val>: non-scan parameter (if flag3=x)"
    exit
end

eps = ARGV[0].to_f

fou = ARGV[1] # which Fourier transform 
[ "c", "s", "p" ].include?( fou ) or
    raise "invalid flag1"

sca = ARGV[2] # scan variable
# slow parameter
slo = case sca
      when "b"
          "w"
      when "w"
          "b"
      else
          raise "invalid flag2"
      end

inp = ARGV[3] # how to specify the slow parameter
# number of scans, fixed value of slow parameter
nsc, val = case inp
           when "r", "s"
               [ ARGV[4].to_i, nil ]
           when "x"
               [ 1, ARGV[4].to_f ]
           else
               raise "invalid flag3"
           end
inp=="s" && nsc<=1 and
    raise "input mode 's' requires at least 2 scans"
inp=="r" and rng = Random.new 

# limits for scan variables and slow parameters
lims = { "b" => [0.1,2.0], "w" => [1e-20,1e10] }
scanlims = lims[sca]
slowlims = lims[slo]

# main loop, with npa values of the slow parameter
(0...nsc).each do |j|
    p = case inp
        when "r"
            x = rng.rand
            exp( (1-x)*log(slowlims[0]) + x*log(slowlims[1]) )
        when "s"
            x = j.to_f/(nsc-1)
            exp( (1-x)*log(slowlims[0]) + x*log(slowlims[1]) )
        when "x"
            val<slowlims[0] || val>slowlims[1] and
                raise "value outside limits"
            val
        end
    scan( eps, fou, sca, slo, p, scanlims )
end
