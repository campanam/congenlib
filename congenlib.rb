#!/usr/bin/env ruby

#-----------------------------------------------------------------------------------------
# congenlib
# Michael G. Campana, 2021
CONGENLIBVER = '0.1.0'
# Smithsonian Conservation Biology Institute
#-----------------------------------------------------------------------------------------

require 'zlib'
require 'bigdecimal'

#-----------------------------------------------------------------------------------------
# Calculate probability density function of chi distribution
def chi_prob(test_val, df)
	if test_val > 0.0
		num = (test_val ** (df/2.0 - 1.0)) * Math.exp(test_val/2.0 * -1.0)
		denom = (2 ** (df/2.0)) * Math.gamma(df/2.0)
		return num/denom
	else
		return 0.0
	end
end
#-----------------------------------------------------------------------------------------
# Determine whether input file is gzipped or not and set method to open it
# From BaitsTools 1.6.8: Campana 2018
def gz_file_open(file) 
	if file[-3..-1] == ".gz"
		yield Zlib::GzipReader.open(file)
	else
		yield File.open(file)
	end
end
#-----------------------------------------------------------------------------------------
# Calculate mean of an array of numbers
# From BaitsTools 1.6.8: Campana 2018
def mean(val = [])
	mean = val.sum.to_f
	mean /= val.size
	return mean
end
#-----------------------------------------------------------------------------------------
# Calculate the standard deviation of an array of numbers
# Modified from CorrSieve 1.7.0: Campana 2011 and FAECES*: Parker et al. 2021
def stdev(val = [], sample = true)
	sample == true ? n = 1 : n = 0 # Determine if sample or population standard deviation
	me = mean(val)
	val.map! { |x| (x.to_f - me) ** 2 }
	st = val.sum
	sd = Math.sqrt( st / (val.size - n))
	return sd
end
#-----------------------------------------------------------------------------------------
# Calculate binomial probability
# From FAECES*: Parker et al. 2021
def binom(n, k, p = 0.50)
	if n - k > k # Maximize efficient division
		nCk = factorial(n, n-k)
		nCk /= factorial(k)
	else
		nCk = factorial(n, k) # Calculate factorial(n)/factorial(k)
		nCk /= factorial(n - k)
	end
	prob = BigDecimal(p.to_s) ** BigDecimal(k.to_s) * BigDecimal((1.0 - p).to_s) ** BigDecimal((n-k).to_s) * BigDecimal(nCk.to_s)
	return prob
end
#-----------------------------------------------------------------------------------------
# Calculate factorials
# From FAECES*: Parker et al. 2021
def factorial(n, k = 1)
	total = 1
	while n > k # Permits efficient division of factorials
		total *= n
		n -= 1
	end
	return total
end
#-----------------------------------------------------------------------------------------
# Calculate chi distribution cumulative probability distribution under special case of df = 1
# From BaitsTools 1.7.3: Campana 2018
def chi_cum_prob(test_val) # This calculates chi distribution cumulative probability distribution under special case of df = 1
	return Math.erf(Math.sqrt(test_val/2.0))
end
#-----------------------------------------------------------------------------------------
# Convert values to BigDecimal for tests
def to_bigdec(vals = [], bigdeclimit = 16)
	BigDecimal.limit(bigdeclimit) # Set Big decimal limit
	return vals.map { |x| BigDecimal(x.to_s) }
end
#-----------------------------------------------------------------------------------------
# Calculate Fisher's exact test
# Modified from amakihi_GO: Paxton et al. in prep
# Gives p value of only exactly this configuration. Not a tailed-test.
def fisher(ref_test, comp_test, ref_total, comp_total, bigdeclimit = 16)
	# ref_test = count of test in the reference
	# comp_test = count of test in the comparison
	# ref_total/comp_total = total values of reference and comparison respectively
	# Convert values to BigDecimal for accurate computation
	ref_test,comp_test,ref_total,comp_total = to_bigdec([ref_test,comp_test,ref_total,comp_total], bigdeclimit)
	not_ref_test = ref_total - ref_test
	not_comp_test = comp_total - comp_test
	denom = factorial(ref_test)*factorial(comp_test)*factorial(not_ref_test)*factorial(not_comp_test)*factorial(ref_total + comp_total)
	numer = factorial(ref_test+comp_test)*factorial(not_ref_test+not_comp_test)*factorial(ref_total)*factorial(comp_total)
	pval = numer/denom
	return pval
end
#-----------------------------------------------------------------------------------------
# Calculate chi2 test
# Modified from amakihi_GO: Paxton et al. in prep
def chi2(ref_test, comp_test, ref_total, comp_total, bigdeclimit = 16)
	# ref_test = count of test in the reference
	# comp_test = count of test in the comparison
	# ref_total/comp_total = total values of reference and comparison respectively
	# Convert values to BigDecimal for accurate computation
	ref_test,comp_test,ref_total,comp_total = to_bigdec([ref_test,comp_test,ref_total,comp_total], bigdeclimit)
	not_ref_test = ref_total - ref_test
	not_comp_test = comp_total - comp_test
	numer = ((ref_test*not_comp_test - comp_test*not_ref_test)**2)*(ref_test+not_comp_test+comp_test+not_ref_test)
	denom = (ref_test+comp_test)*(not_ref_test+not_comp_test)*(comp_test+not_comp_test)*(ref_test+not_ref_test)
	chi2 = numer/denom
	return BigDecimal((1.0 - chi_cum_prob(chi2)).to_s)
end
