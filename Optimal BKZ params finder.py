def H_delta(x):
	#эвристика beta-root Hermite factor
	tmp =( x/(2*pi*e)*(pi*x)^(1/x) )^(1/(2*(x-1)))
	return tmp

def H_delta(x,precision=144):
	#эвристика beta-root Hermite factor
	tmp =( x/(2*pi*e)*(pi*x)^(1/x) )^(1/(2*(x-1))).n(precision)
	return tmp

def left_side(q,n,k):
	return 1/(k*(3*k-1))*(k/2*log(k/2/pi/e,2)+k*log(q,2)-n*log(n/2,2))

q=2^32
n =2^10
def find_param2(n, q, old=False,logs=False):
    min_beta=Infinity
    min_k = int(n/8)
    L=35; R=2*n
    
    while abs(L-R)>1:
        beta=int(L+(R-L)/2)
        
        if logs:
            print("now:",L,beta,R)
            
        delta_beta =log(H_delta(beta),2)
        
        flag_k_found=False
        if not old:
            l=beta; r=2*n

            flag_k_found=False

            while abs(l-r)>1:
                k=int(l+(r-l)/2)
                tmp=left_side(q,n,k)

                #если мы нашли новый минимум k, то поиск продолжается не правее
                if tmp >=delta_beta:
                    flag_k_found=True
                    min_beta = beta
                    min_k = k
                    
                    if logs:
                        print(min_k, min_beta)
                        
                    r=k
                #иначе поиск происходит там, где tmp-delta_beta больше
                else:
                    left_=left_side(q,n,k-1)-delta_beta
                    right_=left_side(q,n,k+1)-delta_beta
                    if left_-right_<0:
                        l=k+1
                    else:
                        r=k
        
        else:
            for k in range(max(beta,min_k), 2*n, 10):
                tmp=left_side(q,n,k)
                if tmp >=delta_beta:
                    min_beta = beta
                    min_k = k
                    print(min_k, min_beta)
                    flag_k_found=True
                    break
        #если для данной бэта есть решение, то ищем новое не правее
        if flag_k_found:
            R=beta
        #иначе ищем решение правее
        else:
            L=beta+1
       
    return min_k, min_beta

print(find_param2(n,q,logs=True))

