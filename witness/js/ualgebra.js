// The MIT License

// Copyright (c) 2011 Pedro http://lamehacks.net

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.



function matrixMultiply(a,b){
	var btrans = transposeMatrix(b);
	var result = [];
	for(var i=0; i < a.length; i++){
		var row = [];
		for(var j=0; j < btrans.length; j++){
			var value = internalProduct(a[i],btrans[j]);
			row.push(value)
		}
		result.push(row);
	}
	return result;
}

function matrixScalarMultiply(m,s){
	var result = [];
	for(var i=0; i < m.length; i++){
		var row = [];
		for(var j=0; j < m[0].length; j++){
			row.push(s * m[i][j]);
		}
		result.push(row);
	}
	return result;
}


//unnecessary - use dotOp instead
function matrixAdd(a,b){
	var result = [];
	for(var i=0; i < a.length; i++){
		var row = [];
		for(var j=0; j < a[0].length; j++){
			row.push(a[i][j] + b[i][j]);
		}
		result.push(row);
	}
	return result;
}

//for internal usage only
function internalProduct(u,v){
	if (u.length != v.length) throw "SizesDoNotMatch";
	var sum = 0;
	for(var i=0; i < u.length; i++){
		sum += u[i]*v[i];
	}
	return sum;
}

function transposeMatrix(m){
	var t = [];
	for(var i=0; i < m[0].length; i++){
		var row = [];
		for(var j=0; j < m.length; j++){
			row.push(m[j][i]);
		}
		t.push(row);
	}
	return t;
}

function minorMatrix(m, k, l){
	var reduced = [];
	for(var i=0; i < m.length; i++){
		if(i==k) continue;
		var row = [];
		for(var j=0; j < m.length; j++){
			if(j==l) continue;
			row.push(m[i][j])
		}
		reduced.push(row);
	}
	return reduced;
}


/*
Recursive implementation using laplace expansion
http://www.webcitation.org/61AGedZlm
*/
function determinant(m){
	var size = m.length;
	if(size == 1) return m[0][0];
	if(size == 2) return m[0][0] * m[1][1] - m[0][1] * m[1][0];
	var det = 0;
	for(var i=0; i < size; i++){
		var minor = minorMatrix(m,0,i);
		var signal = (i%2 > 0) ? -1 : 1;
		det += signal * m[0][i]* determinant(minor);
	}
	return det;
}


/*
http://en.wikipedia.org/wiki/Cofactor_(linear_algebra)
*/
function cofactor(m, k, l){
	minor = minorMatrix(m, k, l);
	return determinant(minor);
}

/*
http://en.wikipedia.org/wiki/Cofactor_(linear_algebra)#Matrix_of_cofactors
*/
function cofactorMatrix(m){
	var cofactors = [];
	for(var i = 0; i < m.length; i++){
		var row = [];
		for(var j = 0; j < m.length; j++){
			var cofactorval = cofactor(m,i,j)*Math.pow(-1,i+j)
			row.push(cofactorval);
		}
		cofactors.push(row);
	}
	return cofactors;
}


/*
Used the adjoint method
http://www.webcitation.org/61BTRqAoZ
*/
function inverseMatrix(m){
	var det = determinant(m);
	if (det == 0) throw "SingularMatrix";
	var deti = 1 / det;
	var cof = cofactorMatrix(m);
	var adj = transposeMatrix(cof);
	var result = matrixScalarMultiply(adj,deti);
	return result;
}

/*
performs operation element by element between to matrices
*/
function dotOp(func,m,n){
	var result = [];
	for(var i = 0; i < m.length; i++){
		var row = [];
		for(var j = 0; j < m[0].length; j++){
			row.push(func(m[i][j], n[i][j]));
		}
		result.push(row);
	}
	return result;
}


function generateMatrix(nlines, ncols, func){
	var m = [];
	for(var i = 0; i < nlines; i++){
		var row = [];
		for(var j = 0; j < ncols; j++){
			row.push(func(i,j));
		}
		m.push(row);
	}
	return m;
}

function zeros(nlines,ncols){
	return generateMatrix(nlines,ncols, function(i,j){return 0;});
}

function ones(nlines,ncols){
	return generateMatrix(nlines,ncols, function(i,j){return 1;});
}

function identity(size){
	return generateMatrix(size, size, function(i,j){ if(i===j) return 1; return 0; });
}
