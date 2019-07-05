stratify = function( Treatment, Outcome, Matrix, Discretize = TRUE , Synthetic = TRUE , Ordered = TRUE, Markov= FALSE) { 

w=0


weightdisc = function (mat) 

{

    if (nrow(mat) == 1) 

        message("Only one strata matched in data, odds ratio may be meaningless on this case.")

    cc = 0

    if (nrow(mat) == 1) 

        cc = 1

    cola = ncol(mat)

    colb = ncol(mat) - 7

    mat = as.matrix(mat[, (colb:cola)])

    if (cc == 1) 

        mat = t(mat)

    sss = sum(mat[, 5])

    l = matrix(nrow = nrow(mat), ncol = ncol(mat))

    for (i in 1:ncol(l)) {

        l[, i] = as.numeric(mat[, i])

    }

    mat = l

    mat[, 3] = mat[, 3] * mat[, 2]

    mat[, 4] = mat[, 4] * mat[, 2]

    n = 0

    a = 0

    b = 0

    c = 0

    d = 0

    run = 0

    run2 = 0

    E = 0

    V = 0

    O = 0

    v1 = 0

    v2 = 0

    v3 = 0

    v4 = 0

    for (i in 1:nrow(mat)) {

        n = mat[i, 3] + mat[i, 4] + mat[i, 7] + mat[i, 8]

        a = mat[i, 7]

        b = mat[i, 8]

        c = mat[i, 3]

        d = mat[i, 4]

        run = run + (a * d)/n

        run2 = run2 + (b * c)/n

        v1 = v1 + (((a * d)/n) * ((a + d)/n))

        v2 = v2 + (((b * c)/n) * ((b + c)/n))

        v3 = v3 + (((a * d)/n) * ((b + c)/n)) + (((b * c)/n) * 

            ((a + d)/n))

        O = O + a

        E = E + ((mat[i, 7] + mat[i, 8]) * (mat[i, 7] + mat[i, 

            3]))/(mat[i, 3] + mat[i, 4] + mat[i, 7] + mat[i, 

            8])

        V = V + ((mat[i, 7] + mat[i, 8]) * (mat[i, 3] + mat[i, 

            8]) * (mat[i, 7] + mat[i, 3]) * (mat[i, 8] + mat[i, 

            4]))/(mat[i, 3] + mat[i, 4] + mat[i, 7] + mat[i, 

            8])^2 * (mat[i, 3] + mat[i, 4] + mat[i, 7] + mat[i, 

            8] - 1)

    }

    v4 = v1/(2 * run^2) + v2/(2 * run2^2) + v3/(2 * run * run2)

    CHI = (abs(O - E) - 0.5)^2/V

    l = pchisq(CHI, df = 1)

    result = run/run2

    k = matrix(nrow = 5, ncol = 1)

    rownames(k) = c("Odds Ratio Of Impact Of Treatment on Outcome", 

        "Mantel Haenszel P-Value", "No. Of Cases Matched", "95% C.I Upper Bound", 

        "95% C.I Lower Bound")

    k[1, 1] = round(result, 3)

    k[2, 1] = round(l, 4)

    k[3, 1] = sss

    k[4, 1] = exp(log(k[1, 1]) + pnorm(0.975, 0, 1) * sqrt(v4))

    k[5, 1] = exp(log(k[1, 1]) - pnorm(0.975, 0, 1) * sqrt(v4))

    return(k)

}

oddsdisc= function (mat) 
{
    if (nrow(mat) == 1) 
        stop("Only one strata matched in data, cannot display odds before and after startification.")
    cola = ncol(mat)
    colb = ncol(mat) - 7
    t = as.matrix(mat[, (colb:cola)])
    t2 = as.matrix(mat[, -(colb:cola)])
    l = matrix(nrow = nrow(t), ncol = ncol(t))
    for (i in 1:ncol(l)) {
        l[, i] = as.numeric(t[, i])
    }
    t = l
    en = matrix(nrow = 1, ncol = 1)
    for (kk in 1:ncol(t2)) {
        d = na.omit(as.numeric(unique(t2[, kk])))
        if (length(d) == 0) {
            s = unique(t2[, kk])
            s = as.matrix(s, byrow = T)
            en = rbind(en, s)
        }
        else {
            s = colnames(t2)[kk]
            a = as.matrix(s)
            en = rbind(en, s)
        }
    }
    en = na.omit(en)
    results = matrix(nrow = nrow(en), ncol = 2)
    results[, 2] = rep(1)
    colnames(results) = c("Before Weighting", "After Weighting")
    count = 1
    for (kk in 1:ncol(t2)) {
        d = na.omit(as.numeric(unique(t2[, kk])))
        if (length(d) == 0) {
            s = unique(t2[, kk])
            for (i in 1:length(s)) {
                dd = as.matrix(t2[, kk])
                for (o in 1:nrow(dd)) {
                  if (dd[o, 1] != s[i]) 
                    dd[o, 1] = 0
                  if (dd[o, 1] == s[i]) 
                    dd[o, 1] = 1
                }
                dd = as.matrix(as.numeric(dd[, 1]))
                a = t[, 5] * dd[, 1]
                b = t[, 1] * dd[, 1]
                a1 = sum(a)/sum(t[, 5])
                b1 = sum(b)/sum(t[, 1])
                f = a1/b1
                results[count, 1] = f
                count = count + 1
            }
        }
        if (length(d) != 0) {
            dd = as.matrix(t2[, kk])
            dd = as.matrix(as.numeric(dd[, 1]))
            a = t[, 5] * dd[, 1]
            b = t[, 1] * dd[, 1]
            a1 = sum(a)/sum(t[, 5])
            b1 = sum(b)/sum(t[, 1])
            f = b1/a1
            results[count, 1] = f
            count = count + 1
        }
    }
    rownames(results) = en[, 1]
    results = na.omit(results)
    return(results)
}
 


dispdisc =function (mat) 
{
    mat = oddsdisc(mat)
    time = c(1:nrow(mat))
    betagal.abs = c(mat[, 1])
    cell.density = c(mat[, 2])
    h = max(mat[, 1])
    h = h + 1
    plot(time, betagal.abs, pch = 17, axes = F, ylim = c(0, h), 
        xlab = "", ylab = "", col = "black", main = "Balance Of Covariates")
    axis(2, ylim = c(0, 1), col = "black")
    mtext("Odds", side = 2, line = 2.5)
    box()
    par(new = T)
    plot(time, cell.density, pch = 15, xlab = "", ylab = "", 
        ylim = c(0, h), axes = F, type = "b", col = "red")
    bb = nrow(mat)
    axis(1, 1:bb, rownames(mat), las = 2)
    legend("topright", legend = c("Before Weighting", "After Weighting"), 
        pch = c(17, 15))
}









weightcont= function (Matrix) 

{

    if (nrow(Matrix) == 1) 

        message("Only one strata matched in data, adjusted regression coefficient may be meaningless in this case.")

    cc = 0

    if (nrow(Matrix) == 1) 

        cc = 1

    resd = matrix(nrow = 5)

    a = ncol(Matrix) - 5

    b = ncol(Matrix)

    r = Matrix[, (a:b)]

    if (cc == 1) {

        r = t(as.matrix(r))

    }

    f = matrix(nrow = nrow(r))

    r = cbind(r, f)

    t = matrix(nrow = nrow(r), ncol = ncol(r))

    for (i in 1:ncol(t)) {

        t[, i] = as.numeric(r[, i])

    }

    r = t

    sum1 = sum(r[, 1])

    r[, 2] = (r[, 1]/sum1) * r[, 2]

    r[, 4] = r[, 1]/r[, 4]

    sum2 = sum(r[, 4])

    r[, 4] = r[, 4]/sum2

    r[, 5] = r[, 5] * r[, 4]

    res = sum(r[, 2]) - sum(r[, 5])

    a = r[, 2] - r[, 5]

    f = try(t.test(a, mu = 0), TRUE)

    if (class(f) == "try-error") {

        resd[4, 1] = NA

        resd[5, 1] = NA

    }

    else {

        resd[4, 1] = f$conf.int[2] * nrow(Matrix)

        resd[5, 1] = f$conf.int[1] * nrow(Matrix)

    }

    resd[1, 1] = res

    resd[2, 1] = as.numeric(f[3])

    resd[3, 1] = round(sum(r[, 1]), 1)

    rownames(resd) = c("Average Of Cases - Average Of Controls", 

        "t-test p-value", "No. Of Cases Matched", "95% C.I Upper Bound", 

        "95% C.I Lower Bound")

    return(resd)

}


oddscont = function (mat) 
{
    if (nrow(mat) == 1) 
        stop("Only one strata matched in data, cannot display odds before and after startification.")
    cola = ncol(mat)
    colb = ncol(mat) - 5
    t = as.matrix(mat[, (colb:cola)])
    t2 = as.matrix(mat[, -(colb:cola)])
    l = matrix(nrow = nrow(t), ncol = ncol(t))
    for (i in 1:ncol(l)) {
        l[, i] = as.numeric(t[, i])
    }
    t = l
    en = matrix(nrow = 1, ncol = 1)
    for (kk in 1:ncol(t2)) {
        d = na.omit(as.numeric(unique(t2[, kk])))
        if (length(d) == 0) {
            s = unique(t2[, kk])
            s = as.matrix(s, byrow = T)
            en = rbind(en, s)
        }
        else {
            s = colnames(t2)[kk]
            a = as.matrix(s)
            en = rbind(en, s)
        }
    }
    en = na.omit(en)
    results = matrix(nrow = nrow(en), ncol = 2)
    results[, 2] = rep(1)
    colnames(results) = c("Before Weighting", "After Weighting")
    count = 1
    for (kk in 1:ncol(t2)) {
        d = na.omit(as.numeric(unique(t2[, kk])))
        if (length(d) == 0) {
            s = unique(t2[, kk])
            for (i in 1:length(s)) {
                dd = as.matrix(t2[, kk])
                for (o in 1:nrow(dd)) {
                  if (dd[o, 1] != s[i]) 
                    dd[o, 1] = 0
                  if (dd[o, 1] == s[i]) 
                    dd[o, 1] = 1
                }
                dd = as.matrix(as.numeric(dd[, 1]))
                a = t[, 4] * dd[, 1]
                b = t[, 1] * dd[, 1]
                a1 = sum(a)/sum(t[, 4])
                b1 = sum(b)/sum(t[, 1])
                f = a1/b1
                results[count, 1] = f
                count = count + 1
            }
        }
        if (length(d) != 0) {
            dd = as.matrix(t2[, kk])
            dd = as.matrix(as.numeric(dd[, 1]))
            a = t[, 4] * dd[, 1]
            b = t[, 1] * dd[, 1]
            a1 = sum(a)/sum(t[, 4])
            b1 = sum(b)/sum(t[, 1])
            f = a1/b1
            results[count, 1] = f
            count = count + 1
        }
    }
    rownames(results) = en[, 1]
    results = na.omit(results)
    return(results)
}



dispcont= function (mat) 
{
    mat = oddscont(mat)
    time = c(1:nrow(mat))
    betagal.abs = c(mat[, 1])
    cell.density = c(mat[, 2])
    h = max(mat[, 1])
    h = h + 1
    plot(time, betagal.abs, pch = 17, axes = F, ylim = c(0, h), 
        xlab = "", ylab = "", col = "black", main = "Balance Of Covariates")
    axis(2, ylim = c(0, 1), col = "black")
    mtext("Odds", side = 2, line = 2.5)
    box()
    par(new = T)
    plot(time, cell.density, pch = 15, xlab = "", ylab = "", 
        ylim = c(0, h), axes = F, type = "b", col = "red")
    bb = nrow(mat)
    axis(1, 1:bb, rownames(mat), las = 2)
    legend("topright", legend = c("Before Weighting", "After Weighting"), 
        pch = c(17, 15))
}



##check data is numeric 
for(i in 1:ncol(Matrix)){
if ( class (Matrix[,i]) != "numeric")stop ("All variables need to be numeric, use as.numeric() function to first change all variables to numeric")
}



##zop is just a counter variable
zop=0

##Make sure the Treatment variable has values either 0 or 1, if not return error 

a=as.numeric(as.character(unique(Matrix[,Treatment])))

if(length(a) != 2) stop ("Treatment variable has problems, either it does not have 2 levels, or it does not take values 0 or 1")

aa=which(a == 0) 

if (length(aa) == 0) stop ("Treatment variable has problems, either it does not have 2 levels, or it does not take values 0 or 1")

aa=which(a == 1) 

if (length(aa) == 0) stop ("Treatment variable has problems, either it does not have 2 levels, or it does not take values 0 or 1")

 

 

##Check if outcome is discrete 

a=as.numeric(as.character(unique(Matrix[,Outcome])))

if(length(a) <= 2){

 

##Make sure the Outcome variable has values either 0 or 1, if not return error 

aa=which(a == 0) 

if (length(aa) == 0) stop ("Outcome variable has problems, either it does not have 2 levels, or it does not take values 0 or 1")

aa=which(a == 1) 

if (length(aa) == 0) stop ("Outcome variable has problems, either it does not have 2 levels, or it does not take values 0 or 1")

 

 

 

##remove all NA's and make sure there is enough data left after cleaning

##if not return error  

Matrix=na.omit(Matrix)

if( nrow(Matrix) < 10) stop("Too few data remaining after removing NA's")

 

 

##check that there are covariates in data, if not return error

if( ncol(Matrix) < 3) stop("No covariates detected in data")

 

##Make sure all the covariates have unique names, if not return error

a=unique(colnames(Matrix))

if(length(a) < ncol(Matrix) ) stop("All variables must have unique names") 

 

 

##orgainze the Matrix so that Treatment and Outcome are the last 2 columns 

mm=cbind(Matrix[,Treatment],Matrix[,Outcome])

colnames(mm) = c( colnames(Matrix)[Treatment], colnames(Matrix)[Outcome])

for(i in 1:ncol(Matrix)){

if(i != Treatment & i != Outcome) {

mm= cbind(Matrix[,i],mm)

colnames(mm)[1]=colnames(Matrix)[i]

}

}

Matrix = mm 

Outcome=ncol(Matrix)

Treatment=ncol(Matrix)-1

 

 

 

##Make sure none of the variables are called "V1" so function ddply 

## will work properly 

a=which(colnames(Matrix) == "V1" ) 

if (length(a) != 0 ) colnames(Matrix)[a] = "Var1" 

 

 

##Convert all values to numeric

for(i in 1:ncol(Matrix) ) { 

Matrix[,i] = as.numeric(as.character((Matrix[,i])))

}

Matrix=as.data.frame(Matrix)


##Discretize the covariates if Discretize = TRUE 

if (Discretize == TRUE ) { 

for(i in 1:ncol(Matrix) ) { 
if(i != Treatment ) {
if(i != Outcome) { 
a=as.numeric(as.character(unique(Matrix[,i])))
if ( length(a) > 10) { 
ss=mean(Matrix[,i])
a=which(Matrix[,i] > ss)
b=which(Matrix[,i] <= ss)
Matrix[a,i]=1
Matrix[b,i]=0
}}}}}


##Convert matrix to data.frame

Matrix=as.data.frame(Matrix)

##Find Markov blanket of treatment if Markov= TRUE 
if ( Markov == TRUE ) {

if ( Ordered == FALSE ) {
print("Since data is unordered, parents of Treatment will be found using Grow-Shrink algorithm")
formula= sprintf("t= gs(Matrix)" )
eval(parse(text=formula))
a=(t$node[[Treatment]]$parents)
z=which(a == colnames(Matrix)[Outcome])
if(length(z) != 0 )a=a[-z]
if(length(a) <= 1) {
print("No parents found using Grow-Shrink, regression will be used instead")
z=glm(Matrix[,Treatment] ~ . , data= Matrix[,(1:(Treatment-1))],family=binomial)
z=as.matrix(summary(z)$coef)
z=z[-1,]
b=which(z[,4] <= 0.05) 
if ( length(b) <= 1 ) {
print("No parents found using regression, only Treatment will be used in stratification")  
zop=1
}
if(length(b) > 1 ){
u=as.matrix(Matrix[,b])
colnames(u) = colnames(Matrix)[b] 
cc=sprintf("Parents found using regression are %s", toString(colnames(u)) )
print(cc) 
Matrix= cbind(u, Matrix[,(Treatment:Outcome)] ) 
}
}
if ( length(a) > 1 ) {
u=as.matrix(Matrix[,a])
cc=sprintf("Parents found using Grow-Shrink are %s", toString(colnames(u)) )
print(cc) 
Matrix= cbind(u, Matrix[,(Treatment:Outcome)] ) 
}}




if ( Ordered == TRUE ) {
print("Since data is ordered, parents of Treatment will be found using regression")
z=glm(Matrix[,Treatment] ~ . , data= Matrix[,(1:(Treatment-1))],family=binomial)
z=as.matrix(summary(z)$coef)
z=z[-1,]
b=which(z[,4] <= 0.05) 
if ( length(b) <= 1 ) {
print("No parents found using regression, only Treatment will be used in stratification")  
zop=1
}
if(length(b) > 1 ){
u=as.matrix(Matrix[,b])
colnames(u) = colnames(Matrix)[b] 
cc=sprintf("Parents found using regression are %s", toString(colnames(u)) )
print(cc) 
Matrix= cbind(u, Matrix[,(Treatment:Outcome)] ) 
}
}


}





##Convert matrix to data.frame

Matrix=as.data.frame(Matrix)
Treatment=ncol(Matrix)-1
Outcome=ncol(Matrix)


if(zop == 1 ) {
oo=glm(Matrix[,Outcome] ~ Matrix[,Treatment],family=binomial)
mm=matrix(nrow=2,ncol=1)
rownames(mm)= c("Odds Ratio Of Impact Of Treatment on Outcome " , "P-Value" )
mm[1,1]=round(exp(oo$coef[2]),3)
mm[2,1]=round(summary(oo)$coef[,4][2],3)
colnames(mm)[1] = "Results"
print(mm)
return(mm)
}

if(zop == 0 ) { 

## FINF UNIQUE STRATA  

a=ncol(Matrix) 

mm=Matrix[,(1:a)]

strata1=ddply(mm, (1:ncol(mm)), nrow)

a=ncol(Matrix) - 1

a1=which(strata1[,a] == 1) 

a2=which(strata1[,a] == 0) 

treat=strata1[a1,]

control=strata1[a2,]

treat2=cbind(treat,treat[,1])

treat2[,ncol(treat2)] = 0

colnames(treat2)[ncol(treat2)] = "Y=1 Controls"

colnames(treat2)[ncol(treat2)-1] = "Y=1 Treatments"

 

##this part matches across treatment variable 

for(kk in 1:nrow(treat)){
w=0
formula= "w= which( "

for( i in 1:ncol(treat)) { 

if(i != Treatment) {

if(i != ncol(treat) ){

if(i != Outcome ) formula = sprintf ( "%s control[,%s] == treat[%s,%s] &" , formula ,i , kk, i) 

if(i == Outcome ) formula = sprintf ( "%s control[,%s] == treat[%s,%s] )" , formula ,i , kk, i) 

}}}

eval(parse(text=formula))

if(length(w) != 0 ) treat2[kk,ncol(treat2)]=control[w,ncol(control)]

if(length(w) == 0 ) treat2[kk,ncol(treat2)]=NA

if(length(w) != 0 ) control[w,ncol(control)] = NA

control=na.omit(control)

}


treat2=na.omit(treat2)

if(nrow(treat2) == 0 ) stop ("No matches found" ) 


##this part matches across outcome variable 

a=which(treat2[,Outcome] == 0) 

treat3=treat2[a,]

treat2=treat2[-a,]  

treat4=cbind(treat2,treat2[,1],treat2[,1])

colnames(treat4)[ncol(treat4)] = "Y=0 Controls"

colnames(treat4)[ncol(treat4)-1] = "Y=0 Treatments"

for(kk in 1:nrow(treat2)){

formula= "w= which( "

for( i in 1:ncol(treat)) { 

if(i < Treatment) {

if(i != (Treatment-1) ) formula = sprintf ( "%s treat3[,%s] == treat2[%s,%s] &" , formula ,i , kk, i) 

if(i == (Treatment-1) ) formula = sprintf ( "%s treat3[,%s] == treat2[%s,%s] )" , formula ,i , kk, i) 

}}

eval(parse(text=formula))

if(length(w) == 0 ) treat4[kk,ncol(treat2)]=NA

if(length(w) != 0 ) treat4[kk,ncol(treat4)]=treat3[w,ncol(treat3)]

if(length(w) != 0 ) treat4[kk,(ncol(treat4)-1)]=treat3[w,(ncol(treat3)-1)]

if(length(w) != 0 ) treat3[w,ncol(treat3)] = NA

treat3=na.omit(treat3)

}



a=which(is.na(treat4[,ncol(treat4) - 2]) == TRUE )  
treat4[a,ncol(treat4) - 2] = 0



treat4=na.omit(treat4)

if(nrow(treat4) == 0 ) stop ("No matches found" ) 

 

treat4=treat4[,-(Treatment:Outcome)]

 

a=ncol(treat4)-3

c=ncol(treat4)-2

b=ncol(treat4) - 1 

d=ncol(treat4)

 

zz=matrix(ncol=8, nrow=nrow(treat4))

colnames(zz) = c('T=0', 'Wk0', 'Y=1', 'Y=0', 'T=1', 'Wk1', 'Y=1' ,'Y=0')

zz[,3]=treat4[,c]

zz[,4]=treat4[,d]

zz[,7]=treat4[,a]

zz[,8]=treat4[,b]

zz[,1]=zz[,3]+zz[,4]

zz[,5]=zz[,7] + zz[,8] 

zz[,6]=rep(1)

zz[,2]= zz[,5]/zz[,1]

infs=which(zz[,2] == Inf)
zz[infs,2]=1

treat4=treat4[,-(a:d)]

treat4=cbind(treat4,zz)

 

 

### this is the function to do the weighting 


 

o=weightdisc(treat4)
castrino=o[3,1]/ length(which(Matrix[,Treatment] == 1 )) 
if(Synthetic == TRUE ) {
l=glm(Matrix[,ncol(Matrix)] ~ . ,data = Matrix[,(1:(ncol(Matrix)-1))], family=binomial)
a=which(Matrix[,Treatment] == 1) 
aa=l$coef
aa=aa[length(aa)]
weightSCB=o[3,1]/length(a)
weightGLM=1-weightSCB
sss=weightSCB*o[1,1] + weightGLM*exp(aa)
o[1,1]=sss
o[3,1]=as.numeric(length(a))
castrino2=o[3,1]/ length(which(Matrix[,Treatment] == 1 )) 
zz=sprintf ("Before usning synthetic mathcing, percent of cases mathced was: %s", castrino*100) 
zz2=sprintf ("After usning synthetic mathcing, percent of cases mathced became: %s", castrino2*100) 
print(zz)
print(zz2)
}

for(i in 1:nrow(o)) {
o[i,1] = round(o[i,1] , 3)
} 
if(o[2,1] > 0.5) o[2,1] = 1- 0.5 
o[,1]=round(o[,1],3)


dispdisc(treat4)

 

ll=list()

ll[[1]]=treat4
names(ll)[[1]]="Strata"

ll[[2]]=o
names(ll)[[2]]="OddsRatio"

colnames(o)[1] = "Results"
print(o)

return(ll)

}


if( zop ==1 ) print(mm)
}

 

 

## CHECK IF OUTCOME IS CONTINUOUS 

a=as.numeric(as.character(unique(Matrix[,Outcome])))

if(length(a) > 2){

 

##remove all NA's and make sure there is enough data left after cleaning

##if not return error  

Matrix=na.omit(Matrix)

if( nrow(Matrix) < 10) stop("Too few data remaining after removing NA's")

 

 

##check that there are covariates in data, if not return error

if( ncol(Matrix) < 3) stop("No covariates detected in data")

 

##Make sure all the covariates have unique names, if not return error

a=unique(colnames(Matrix))

if(length(a) < ncol(Matrix) ) stop("All variables must have unique names") 

 

 

##orgainze the Matrix so that Treatment and Outcome are the last 2 columns 

mm=cbind(Matrix[,Treatment],Matrix[,Outcome])

colnames(mm) = c( colnames(Matrix)[Treatment], colnames(Matrix)[Outcome])

for(i in 1:ncol(Matrix)){

if(i != Treatment & i != Outcome) {

mm= cbind(Matrix[,i],mm)

colnames(mm)[1]=colnames(Matrix)[i]

}

}

Matrix = mm 

Outcome=ncol(Matrix)

Treatment=ncol(Matrix)-1

 

 

 

##Make sure none of the variables are called "V1" so function ddply 

## will work properly 

a=which(colnames(Matrix) == "V1" ) 

if (length(a) != 0 ) colnames(Matrix)[a] = "Var1" 

 

 

##Convert all values to numeric

for(i in 1:ncol(Matrix) ) { 

Matrix[,i] = as.numeric(as.character((Matrix[,i])))

}

 

##Convert matrix to data.frame

Matrix=as.data.frame(Matrix)



##Convert matrix to data.frame

Matrix=as.data.frame(Matrix)


##Discretize the covariates if Discretize = TRUE 

if (Discretize == TRUE ) { 

for(i in 1:ncol(Matrix) ) { 
if(i != Treatment ) {
if(i != Outcome) { 
a=as.numeric(as.character(unique(Matrix[,i])))
if ( length(a) > 10) { 
ss=mean(Matrix[,i])
a=which(Matrix[,i] > ss)
b=which(Matrix[,i] <= ss)
Matrix[a,i]=1
Matrix[b,i]=0
}}}}}






##Find Markov blanket of treatment if Markov= TRUE 


if ( Markov == TRUE ) {

if ( Ordered == FALSE ) {
print("Since data is unordered, parents of Treatment will be found using Grow-Shrink algorithm")
formula= sprintf("t= gs(Matrix)" )
eval(parse(text=formula))
a=(t$node[[Treatment]]$parents)
z=which(a == colnames(Matrix)[Outcome])
if(length(z) != 0 )a=a[-z]
if(length(a) <= 1) {
print("No parents found using Grow-Shrink, regression will be used instead")
z=glm(Matrix[,Treatment] ~ . , data= Matrix[,(1:(Treatment-1))],family=binomial)
z=as.matrix(summary(z)$coef)
z=z[-1,]
b=which(z[,4] <= 0.05) 
if ( length(b) <= 1 ) {
print("No parents found using regression, all variables will be used in stratification")  
zop=1
}
if(length(b) > 1 ){
u=as.matrix(Matrix[,b])
colnames(u) = colnames(Matrix)[b] 
cc=sprintf("Parents found using regression are %s", toString(colnames(u)) )
print(cc) 
Matrix= cbind(u, Matrix[,(Treatment:Outcome)] ) 
}
}
if ( length(a) > 1 ) {
u=as.matrix(Matrix[,a])
colnames(u) = colnames(Matrix)[b] 
cc=sprintf("Parents found using Grow-Shrink are %s", toString(colnames(u)) )
print(cc) 
Matrix= cbind(u, Matrix[,(Treatment:Outcome)] ) 
}}




if ( Ordered == TRUE ) {
print("Since data is ordered, parents of Treatment will be found using regression")
z=glm(Matrix[,Treatment] ~ . , data= Matrix[,(1:(Treatment-1))],family=binomial)
z=as.matrix(summary(z)$coef)
z=z[-1,]
b=which(z[,4] <= 0.05) 
if ( length(b) <= 0 ){
print("No parents found using regression, all variables will be used in stratification")  
zop=1
}
if(length(b) > 1 ){
u=as.matrix(Matrix[,b])
colnames(u) = colnames(Matrix)[b] 
cc=sprintf("Parents found using regression are %s", toString(colnames(u)) )
print(cc) 
Matrix= cbind(u, Matrix[,(Treatment:Outcome)] ) 
}
}


}






##Convert matrix to data.frame

Matrix=as.data.frame(Matrix)

Treatment=ncol(Matrix)-1
Outcome=ncol(Matrix)


if(zop == 1 ) {
a=which(Matrix[,Treatment] == 1) 
b=which(Matrix[,Treatment] == 0)
a1=mean(Matrix[a,Outcome])
b1=mean(Matrix[b,Outcome])
mm=matrix(nrow=2,ncol=1)
rownames(mm)= c("Average of Cases - Avverage of Controls" , "P-Value" )
mm[1,1]=round(a1-b1,3)
sss=try(t.test(Matrix[a,Treatment],Matrix[b,Treatment])$p.value,silent =TRUE)
if(class(sss) == "try-error")mm[2,1]=1
if(class(sss) != "try-error")mm[2,1]=sss
colnames(mm)[1] = "Results"
print(mm)
return(mm)
}

if(zop ==0) {
## FIND UNIQUE STRATA  

a=ncol(Matrix) - 1 

mm=Matrix[,(1:a)]

strata1=ddply(mm, (1:ncol(mm)), nrow)

a=ncol(Matrix) - 1

a1=which(strata1[,a] == 1) 

a2=which(strata1[,a] == 0) 

treat=strata1[a1,]

control=strata1[a2,]

treat=cbind(treat,treat[,1])

colnames(treat)[ncol(treat)] = "AverageTreatment"

control=cbind(control,control[,1])

colnames(control)[ncol(control)] = "AverageControl"

 

##FIND THE AVERAGES FOR TREATMENTS

for(kk in 1:nrow(treat)){

formula= "w= which( "

for( i in 1:ncol(treat)) { 

if(i != Outcome) {

if(i != ncol(treat) ){

if(i != Treatment ) formula = sprintf ( "%s Matrix[,%s] == treat[%s,%s] &" , formula ,i , kk, i) 

if(i == Treatment ) formula = sprintf ( "%s Matrix[,%s] == treat[%s,%s] )" , formula ,i , kk, i) 

}}}

eval(parse(text=formula))

if(length(w) != 0 ) treat[kk,ncol(treat)]=mean(Matrix[w,ncol(Matrix)]) 

if(length(w) == 0 ) treat[kk,ncol(treat)]=NA

}

treat=na.omit(treat)

if(nrow(treat) == 0 ) stop ("No matches found" ) 

 

##FIND THE AVERAGES FOR CONTROLS

for(kk in 1:nrow(control)){

formula= "w= which( "

for( i in 1:ncol(control)) { 

if(i != Outcome) {

if(i != ncol(control) ){

if(i != Treatment ) formula = sprintf ( "%s Matrix[,%s] == control[%s,%s] &" , formula ,i , kk, i) 

if(i == Treatment ) formula = sprintf ( "%s Matrix[,%s] == control[%s,%s] )" , formula ,i , kk, i) 

}}}

eval(parse(text=formula))

if(length(w) != 0 ) control[kk,ncol(control)]=mean(Matrix[w,ncol(Matrix)]) 

if(length(w) == 0 ) control[kk,ncol(control)]=NA

}

control=na.omit(control)

if(nrow(control) == 0 ) stop ("No matches found" )

 

 

gg=matrix(ncol=4,nrow=nrow(treat),data=0) 

#MATCH TREATMENT AND CONTROL

bb=Treatment -1  

for(kk in 1:nrow(treat)){

formula= "w= which( "

for( i in 1:bb) { 

if(i != ncol(treat) ){

if(i != bb ) formula = sprintf ( "%s control[,%s] == treat[%s,%s] &" , formula ,i , kk, i) 

if(i == bb ) formula = sprintf ( "%s control[,%s] == treat[%s,%s] )" , formula ,i , kk, i) 

}}

eval(parse(text=formula))

if(length(w) != 0 ) {

gg[kk,2]=control[w,ncol(control)-1]

gg[kk,3]=control[w,ncol(control)]

}

if(length(w) == 0 ) gg[kk,4]=NA

}

treat=cbind(treat,gg)

treat=na.omit(treat)

if(nrow(treat) == 0 ) stop ("No matches found" ) 

 

## this is the function to do the weighting 

## do the weighting 

o=weightcont(treat)
dispcont(treat)
 

 

ll=list()

treat=treat[,-ncol(treat)]
treat=treat[,-(ncol(treat)-2)]
colnames(treat)[ncol(treat)-3]='No.OfCases' 
colnames(treat)[ncol(treat)-1]='No.OfControls' 
colnames(treat)[ncol(treat)]='AverageControls'
colnames(treat)[ncol(treat)-2]='AverageCases'




ll[[1]]=treat

names(ll)[[1]]="Strata"

if(o[2,1] > 0.5) o[2,1] = 1- 0.5 
o[,1]=round(o[,1],3)

ll[[2]]=o

names(ll)[[2]]="AverageTreatmentEffect"

colnames(o)[1] = "Results"

print(o)

return(ll)

}}

 

## this parenthesis closes the big function stratify() 

}






