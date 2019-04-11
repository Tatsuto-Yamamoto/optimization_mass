l = ["up", "charm","down","strange","botom","electron","muon","tauon"]
b = ""
for i in l:
	# b+='-2 * (1 - pmas.{}s(x)) * pmas.d{}s(x)'.format(i,i)
	# b += 'np.log(pmas.{}s(x))**2+'.format(i,i)
	b += '2*np.log(pmas.{}s(x))*pmas.d{}s(x)/pmas.{}s(x)+'.format(i,i,i)
print(b)