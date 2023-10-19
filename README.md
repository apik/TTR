# TTR - Tadpoles Tensor Reduction


1. In the input expression all indices should be unique

2. Produced result contains all indices `mu` replaced with `mu?` for FORM substitution

3. All rational coefficients are wrapped with function `rat(a,b) = a/b`

### Usage:

Print results to standard output

```
$ echo "p1(mu1)*p1(mu2)*p1(mu3)*p1(mu4)*p2(mu5)*p2(mu6)" | ./TTR
```


Save to the file FILENAME

```
$ echo "p1(mu1)*p1(mu2)*p1(mu3)*p1(mu4)*p2(mu5)*p2(mu6)" | ./TTR -o FILENAME
```

### Example session

Input for the compact output production:

```
echo "p1(mu1)*p2(mu2)*p3(mu3)*p4(mu4)" | ./TTR -c -o ex4r
```

Content of the file `ex4r`:

```
* 
* 
* [p1.mu1 p2.mu2 p3.mu3 p4.mu4]
* 
* 
*--#[ declare :
S ttrB1;
S ttrB2;
*--#] declare :
*--#[ reduce :
 id p1(mu1?)*p2(mu2?)*p3(mu3?)*p4(mu4?) = 

 + d_(mu1,mu3)*d_(mu2,mu4) * (
	+ p1.p2*p3.p4 * ttrB2 * 1
	+ p1.p4*p2.p3 * ttrB2 * 1
	+ p1.p3*p2.p4 * ttrB1 * 1
   )

 + d_(mu1,mu4)*d_(mu2,mu3) * (
	+ p1.p2*p3.p4 * ttrB2 * 1
	+ p1.p4*p2.p3 * ttrB1 * 1
	+ p1.p3*p2.p4 * ttrB2 * 1
   )

 + d_(mu1,mu2)*d_(mu3,mu4) * (
	+ p1.p2*p3.p4 * ttrB1 * 1
	+ p1.p4*p2.p3 * ttrB2 * 1
	+ p1.p3*p2.p4 * ttrB2 * 1
   )

;
*--#] reduce :
*--#[ subs :
 id ttrB1 = rat(+ (-1/2)+ (-1/2)*d^1,+(-1/2)*d^3+(-1/2)*d^2+(1/1)*d^1);
 id ttrB2 = rat(+ (1/1),+(-1/2)*d^3+(-1/2)*d^2+(1/1)*d^1)/2;
*--#] subs :
```
