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
## This script determines curves omega(beta) that delimit the domains
## where kwws and kwwc can be computed using the low- or high-w expansion.
##
##############################################################################


#! Call runkww with given arguments, check for errors, return kwws or kwwc.

def callkww( argstr )
    cmd = "runkww 0 #{argstr} > /ram/b 2> /ram/c"
    system( cmd )
    c = `cat /ram/c`
    c != "" and
        raise "command '#{cmd}' yields invalid result:\n#{c}"
    return `cat /ram/b`.split(/\n/)[0].to_f
end


#! Determine one limiting curve.

#  input: 
#    dir   = which transform
#    alg   = which series expansion
#    start = w-value for which the chosen algorithm certainly works
#  returns:
#    an array with ln(omega) values

def findlim( dir, alg, start )
    ret = []      # return array
    fac = 2.0     # multiplicative step for initial search
    up = alg=='l' # search with ascending omega
    $bvec.each do |b|
        # 'start' must be in the domain of 'alg'
        wfail = start
        r = callkww( "#{dir} #{alg} #{b} #{wfail}" )
        r >= 0 or
            raise "invalid starting value for b=#{b}: #{wfail} -> '#{r}'"
        # initial search for w outside the domain of 'alg'
        while true
            wgood = wfail
            wfail *= ( up ? fac : 1/fac )
            ( wfail<1e-40 || wfail>1e40 ) and break
            r = callkww( "#{dir} #{alg} #{b} #{wfail}" )
            r < 0 and break
        end
        ( wfail<1e-40 || wfail>1e40 ) and
            raise "#{dir} #{b} no limit"
        # bisection on logarithmic w scale, to find limit of 'alg' domain
        wmid = nil
        30.times do
            wmid = exp( ( log(wfail) + log(wgood) ) / 2 )
            r = callkww( "#{dir} #{alg} #{b} #{wmid}" )
            failed = r < 0
            if up
                wgood,wfail = failed ? [ wgood,wmid ] : [ wmid,wfail ]
            else
                wfail,wgood = failed ? [ wmid,wgood ] : [ wfail,wmid ]
            end
        end
        # value found: save in return array, and print progress report
        ret.push( wmid )
        puts "#{b} #{wmid}"
    end
    return ret
end


#! Main program.

unless ARGV.size==1
    puts "usage:"
    puts "   kww_findlimits.rb c|s|p"
    exit
end

dir = ARGV[0]

# set number of beta values for which limits will be determined
$nbeta = 800

# set logarithmic grid of beta values
$bvec = (0..$nbeta).collect{ |i| 0.1 * (1.999/0.1)**(i.to_f/$nbeta) }

# make sure output files are writable before starting the computation
ftab = File.open( "limits-#{dir}.tab", "w" )

# the computation
puts "low-w limit"
liml = findlim( dir, 'l', 1.0e-20 )
puts "high-w limit"
limh = findlim( dir, 'h', 1.0e+10 )

# write tabular output
(0..$nbeta).each do |i|
    ftab.printf( "%15.9g %15.9g %15.9g\n" % [ $bvec[i], liml[i], limh[i] ] )
end
ftab.close
