# Bond length statistics from X-ray diffraction data
# S. Hessam M. Mehr, 2017, MIT license

# Terrible terrible code that happens to work

module XRD

using Cubature:hcubature

gaussian(mu,sigma) = x-> exp(-(x-mu)^2/(2sigma^2))/sqrt(2pi*sigma^2)

normalize_p(f,min,max) = begin a = hquadrature(f,min,max)[1]; x->f(x)/a end

mean_p(f,min,max) = begin ff = normalize_p(f,min,max); hquadrature(x->x*ff(x),min,max)[1] end

std_p(f,min,max) = begin ff=normalize_p(f,min,max); μ = mean_p(ff,min,max); sqrt(hquadrature(x->(x-μ)^2*ff(x),min,max)[1]) end

function parse_bond(line)
    v1,v2 = split(strip(line),"(")
    v3=length(split(v1,".")[2])
    (float(v1),float(v2[1:end-1])*10.0^-v3) 
end

function distribution(line)
    b=parse_bond(line)
    gaussian(b[1],b[2])
end

function stats(lines, min, max; delim=",") 
    l = split(lines,delim)
    ds = map(distribution,l)
    overall = x->sum(d(x) for d in ds)
    (mean_p(overall,min,max),std_p(overall,min,max))
end

end
