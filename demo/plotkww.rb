#!/usr/bin/env ruby

unless ARGV.size==2
    puts "usage:"
    puts "   plotkww c|s|p l|m|h"
    exit
end

bvec = [ 0.3, 0.45, 0.65, 1, 1.5, 2 ]
 
nw = 500
wl = 1e-5
wh = 1e3
wvec = (0..nw).collect{ |i| wl * (wh/wl)**(i.to_f/(nw-1)) }

bvec.each do |b|
    puts "#{b}"
    wvec.each do |w|
        r = `runkww 0 #{ARGV[0]} #{ARGV[1]} #{b} #{w}`
        a = r.split()
        v = a[0].to_f
        v>0 or next
        puts "#{w} #{a[0]}"
    end
    puts ""
end

use_this_for_copy_and_paste = <<EOT
~/code/kww/test/plotkww.rb c m > hcm.tab
~/code/kww/test/plotkww.rb c l > hcl.tab
~/code/kww/test/plotkww.rb c h > hch.tab
~/code/kww/test/plotkww.rb s h > hsh.tab
~/code/kww/test/plotkww.rb s l > hsl.tab
~/code/kww/test/plotkww.rb s m > hsm.tab
~/code/kww/test/plotkww.rb p m > hpm.tab
~/code/kww/test/plotkww.rb p l > hpl.tab
~/code/kww/test/plotkww.rb p h > hph.tab

ftvm hsl.tab 
ftvm hsm.tab 
ftvm hsh.tab 
ftvm hcl.tab 
ftvm hcm.tab 
ftvm hch.tab 
ftvm hpl.tab 
ftvm hpm.tab 
ftvm hph.tab 

cca t/(1+t^2)
cca dawson(t/2)
cca 1/(1+t^2)
cca sqrt(pi)/2*exp(-t^2/4)

EOT
