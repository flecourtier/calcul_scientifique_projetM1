/*
Factorisation d'une matrice tridiagonale
*/

// A matrice tridiagonale avec low sa sous-diagonale, diag sa diagonale et sup sa sur-diagonale

//calcul de la factorisation LU de la matrice
pub fn factolu(mut l: Vec<f64>, mut d: Vec<f64>, u: Vec<f64>) 
-> (Vec<f64>, Vec<f64>, Vec<f64>){
    let n = d.len();
    for i in 0..n - 1 {
        l[i] = l[i] / d[i];
        d[i + 1] = d[i + 1] - l[i] * u[i];
    }
    (l,d,u)
}

//algorithme de descente-remontée pour résoudre Ax=b
//on suppose qu'un appel précédent a permis de factoriser A
pub fn lusolve(l: &Vec<f64>, d: &Vec<f64>, u: &Vec<f64>, mut x: Vec<f64>, b: &Vec<f64>)
-> Vec<f64> {
    x[0] = b[0];
    let n = d.len();
    // descente
    for i in 1..n {
        x[i] = b[i] - l[i - 1] * x[i - 1];
    }
    //remontée
    x[n - 1] = x[n - 1] / d[n - 1];
    for i in (0..n - 1).rev() {
        x[i] = (x[i] - u[i] * x[i + 1]) / d[i];
    }
    x
}

//produit Ax et stocke le résultat dans b
//on suppose qu'un appel précédent a permis de factoriser A
// pub fn produit_matvec(l: &Vec<f64>, d: &Vec<f64>, u: &Vec<f64>, x: &Vec<f64>, mut b: Vec<f64>)
// -> Vec<f64> {
//     let n = d.len();

//     //Pour avoir la matrice A à partir de la facto LU
//     let mut l_a = vec![0.;n-1];
//     let mut d_a = vec![0.;n];
//     let mut u_a = vec![0.;n-1];
//     d_a[0]=d[0];
//     for i in 0..n-2 {
//         l_a[i]=l[i]*d[i]; 
//         d_a[i+1]=l[i]*u[i]+d[i+1];
//         u_a[i]=u[i];
//     }

//     b[0]=d[0]*x[0]+u[0]*x[1];
//     for i in 1..n-2 {
//         b[i] = l_a[i-1]*x[i-1]+d_a[i]*x[i]+u_a[i]*x[i+1];
//     }
//     b[n-1] = l_a[n-2]*x[n-2]+d_a[n-1]*x[n-1];
//     b
// }

pub fn produit_matvec(l: &Vec<f64>, d: &Vec<f64>, u: &Vec<f64>, x: &Vec<f64>, mut b: Vec<f64>)
-> Vec<f64> {
    let n = d.len();

    //Pour avoir la matrice A e partir de la facto LU
    let mut l_a = vec![0.;n-1];
    let mut d_a = vec![0.;n];
    let mut u_a = vec![0.;n-1];
    d_a[0]=d[0];
    for i in 0..n-1 {
        l_a[i]=l[i]*d[i];
        d_a[i+1]=l[i]*u[i]+d[i+1];
        u_a[i]=u[i];
    }

    b[0]=d[0]*x[0]+u[0]*x[1];
    for i in 1..n-1 {
        b[i] = l_a[i-1]*x[i-1]+d_a[i]*x[i]+u_a[i]*x[i+1];
    }
    b[n-1] = l_a[n-2]*x[n-2]+d_a[n-1]*x[n-1];
    b

    // b[0]=d[0]*x[0]+u[0]*x[1];
    // for i in 1..n-1 {
    //     b[i] = l[i-1]*x[i-1]+d[i]*x[i]+u[i]*x[i+1];
    // }
    // b[n-1] = l[n-2]*x[n-2]+d[n-1]*x[n-1];
    // b
}