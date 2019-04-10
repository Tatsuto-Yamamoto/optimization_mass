l = ["up", "charm","down","strange","botom","electron","muon","tauon"]
b = ""
for i in l:
	b+='-2 * (1 - pmas.{}s(x)) * pmas.d{}s(x)'.format(i,i)
print(b)