use std::f64::consts::PI;

//QUETION 3 

fn define_A(L: f64, N: usize, M:usize) -> Vec<(usize,usize,f64)> {
    let h = L/((N+1) as f64); 

    let taille_A=(N+2)*(M+2);

    //on initialise les 3 tableaux vides 
    //vecval = (i,j,val)
    let mut vecval = vec![];

    //on initialise les valeurs de la diagonales
    for k in 0..taille_A {
        let i = k % (N + 2);
        let j = k / (N + 2);
        if i == 0 || i == N || j == 0 || j == M {
            vecval.push((k, k, 1e20));
        }
        else {
            vecval.push((k, k, 4./h/h));
        }
    }

    //on initialise la sous_diagonale et la sure diagonale
    for i in 0..N+1 {
        for j in 0..M+2 {
            let k1 = j * (N + 2) + i;
            let k2 = j * (N + 2) + i + 1;
            vecval.push((k1, k2, -1. / h / h));
            vecval.push((k2, k1, -1. / h / h));
        }
    }

    //on initialise la 3ème sous_diagonale et la 3ème sure diagonale
    for j in 0..M+1 {
        for i in 0..N + 2 {
            let k1 = j * (N + 2) + i;
            let k2 = (j+1) * (N + 2) + i;
            vecval.push((k1, k2, -1. / h / h));
            vecval.push((k2, k1, -1. / h / h));
        }
    }

    vecval
}

// QUESTION 6

fn define_A_plein(L: f64, N: usize, M:usize) -> Vec<(usize,usize,f64)> {
    let h = L/((N+1) as f64); 

    let taille_A=(N+2)*(M+2);

    //on initialise vecval avec un stockage plein 
    //vecval = (i,j,val)
    let mut vecval = vec![];

    //on initialise les valeurs de la diagonales
    for k in 0..taille_A {
        let i = k % (N + 2);
        let j = k / (N + 2);
        if i == 0 || i == N || j == 0 || j == M {
            vecval.push((k, k, 1e20));
        }
        else {
            vecval.push((k, k, 4./h/h));
        }
    }

    //on initialise la sous_diagonale et la sure diagonale
    for i in 0..N+1 {
        for j in 0..M+2 {
            let k1 = j * (N + 2) + i;
            let k2 = j * (N + 2) + i + 1;
            vecval.push((k1, k2, -1. / h / h));
            vecval.push((k2, k1, -1. / h / h));
        }
    }

    //on initialise la 3ème sous_diagonale et la 3ème sure diagonale
    for j in 0..M+1 {
        for i in 0..N + 2 {
            let k1 = j * (N + 2) + i;
            let k2 = (j+1) * (N + 2) + i;
            vecval.push((k1, k2, -1. / h / h));
            vecval.push((k2, k1, -1. / h / h));
        }
    }

    //on crée le stockage plein 
    for i in 0..taille_A {
        for j in 0..taille_A {
            vecval.push((i, j, 0.));
        }
    }

    vecval
}

//QUESTION 5

//on définit la fonction f(x,y)
fn f(x:f64, y:f64, L: f64, N: usize, M:usize) -> f64 {
    let H=L/((N+1) as f64)*((M+1) as f64);
    PI*PI*(1./(L*L)+1./(H*H))*f64::sin(PI/L*x)*f64::sin(PI/H*y)
}

//on définit u(x,y)
fn u(x:f64, y:f64, L: f64, N: usize, M:usize) -> f64 {
    let H=L/((N+1) as f64)*((M+1) as f64);
    f64::sin(PI/L*x)*f64::sin(PI/H*y)
}



// QUESTION 4
fn define_F(L: f64, N: usize, M:usize) -> Vec<f64> {
    let h = L/((N+1) as f64); 
    // let H = h*((M+1) as f64);
    let taille_F=(N+2)*(M+2);
    let mut F =vec![0.;taille_F];
    for k in 0..taille_F {
        let i = k%(N+2);
        let j = k/(N+2);
        F[k]=f((i as f64)*h, (j as f64)*h, L, N, M);
    }
    F
}

fn solve_laplace2D(L: f64, N: usize, M: usize) -> Vec<f64> {
    let mut A = skyrs::Sky::new(define_A_plein(L,N,M));
    let mut F = define_F(L, N, M);
    let sol = A.solve(F).unwrap();
    sol
}

/// Plot a 2D data set using matplotlib
fn plotpy(xp: Vec<f64>, yp: Vec<f64>, zp: Vec<f64>) {
    use std::fs::File;
    use std::io::BufWriter;
    use std::io::Write;
    {
        let meshfile = File::create("plotpy.dat").unwrap();
        let mut meshfile = BufWriter::new(meshfile);
        xp.iter().for_each(|x| {
            writeln!(meshfile, "{}", x).unwrap();
        });
        writeln!(meshfile, "").unwrap();
        yp.iter().for_each(|y| {
            writeln!(meshfile, "{}", y).unwrap();
        });
        writeln!(meshfile, "").unwrap();
        zp.iter().for_each(|z| {
            writeln!(meshfile, "{}", z).unwrap();
        });
    }

    use std::process::Command;
    Command::new("python3")
        .arg("src/plot.py")
        .status()
        .expect("plot failed !");
}


fn main() {
    let L=1.;
    let N=100;
    let M=100;
    let h = L/((N+1) as f64);
    let xp: Vec<f64> = (0..N + 2).map(|i| i as f64 * h).collect();
    let yp: Vec<f64> = (0..M + 2).map(|i| i as f64 * h).collect();
    let sol = solve_laplace2D(L, N, M);
    // plotpy(xp, yp, sol);
}
