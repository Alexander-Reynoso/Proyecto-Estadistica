// ===============================
// FUNCIONES DE APOYO
// ===============================

// Factorial de un número (n!)
function factorial(num) {
    if (num < 0) return NaN;
    if (num === 0 || num === 1) return 1;
    let result = 1;
    for (let i = 2; i <= num; i++) {
        result *= i;
    }
    return result;
}

// ===============================
// PERMUTACIONES
// ===============================
function calcularPermutacion() {
    const n = parseInt(document.getElementById("nPerm").value);
    const r = parseInt(document.getElementById("rPerm").value);

    if (isNaN(n) || isNaN(r) || n < 0 || r < 0 || r > n) {
        document.getElementById("resultadoPerm").innerText = "⚠️ Por favor, ingresa valores válidos (n ≥ r ≥ 0).";
        return;
    }

    const resultado = factorial(n) / factorial(n - r);
    document.getElementById("resultadoPerm").innerText = `Resultado: ${resultado.toLocaleString()}`;
    document.getElementById("interpretacionPerm").innerText =
        `Existen ${resultado.toLocaleString()} formas diferentes de ordenar ${r} elementos tomados de un total de ${n}.`;
}

// ===============================
// COMBINACIONES
// ===============================
function calcularCombinacion() {
    const n = parseInt(document.getElementById("nComb").value);
    const r = parseInt(document.getElementById("rComb").value);

    if (isNaN(n) || isNaN(r) || n < 0 || r < 0 || r > n) {
        document.getElementById("resultadoComb").innerText = "⚠️ Ingresa valores válidos (n ≥ r ≥ 0).";
        return;
    }

    const resultado = factorial(n) / (factorial(r) * factorial(n - r));
    document.getElementById("resultadoComb").innerText = `Resultado: ${resultado.toLocaleString()}`;
    document.getElementById("interpretacionComb").innerText =
        `Existen ${resultado.toLocaleString()} maneras diferentes de seleccionar ${r} elementos de un total de ${n}.`;
}

// ===============================
// PROBABILIDAD SIMPLE
// ===============================
function calcularProbabilidad() {
    const favorables = parseFloat(document.getElementById("favorables").value);
    const posibles = parseFloat(document.getElementById("posibles").value);

    if (isNaN(favorables) || isNaN(posibles) || favorables < 0 || posibles <= 0 || favorables > posibles) {
        document.getElementById("resultadoProb").innerText = "⚠️ Valores inválidos. Asegúrate de que 0 ≤ favorables ≤ posibles.";
        return;
    }

    const prob = favorables / posibles;
    document.getElementById("resultadoProb").innerText = `P(A) = ${prob.toFixed(3)}`;
    document.getElementById("interpretacionProb").innerText =
        `La probabilidad de que ocurra el evento A es de ${(prob * 100).toFixed(1)}%.`;
}

// ===============================
// EVENTOS INDEPENDIENTES
// ===============================
function calcularIndependientes() {
    const pA = parseFloat(document.getElementById("probA").value);
    const pB = parseFloat(document.getElementById("probB").value);

    if (isNaN(pA) || isNaN(pB) || pA < 0 || pB < 0 || pA > 1 || pB > 1) {
        document.getElementById("resultadoInd").innerText = "⚠️ Ingresa probabilidades válidas (entre 0 y 1).";
        return;
    }

    const pAB = pA * pB;
    document.getElementById("resultadoInd").innerText = `P(A ∩ B) = ${pAB.toFixed(3)}`;
    document.getElementById("interpretacionInd").innerText =
        `La probabilidad de que ambos eventos ocurran simultáneamente es ${(pAB * 100).toFixed(1)}%.`;
}

// ===============================
// EVENTOS EXCLUYENTES (REGLA DE LA SUMA)
// ===============================
function calcularSuma() {
    const pA = parseFloat(document.getElementById("probA2").value);
    const pB = parseFloat(document.getElementById("probB2").value);

    if (isNaN(pA) || isNaN(pB) || pA < 0 || pB < 0 || pA > 1 || pB > 1) {
        document.getElementById("resultadoSuma").innerText = "⚠️ Ingresa probabilidades válidas (entre 0 y 1).";
        return;
    }

    let pUnion = pA + pB;
    if (pUnion > 1) pUnion = 1; // Evita pasar el 100%

    document.getElementById("resultadoSuma").innerText = `P(A ∪ B) = ${pUnion.toFixed(3)}`;
    document.getElementById("interpretacionSuma").innerText =
        `La probabilidad de que ocurra A o B (o ambos) es ${(pUnion * 100).toFixed(1)}%.`;
}

// ----- PERMUTACIONES CON REPETICIÓN -----
function calcularPermutacionRepeticion() {
    let n = parseInt(document.getElementById("nTotal").value);
    let repeticionesTexto = document.getElementById("repeticiones").value.trim();

    if (isNaN(n) || n <= 0 || repeticionesTexto === "") {
        document.getElementById("resultadoPermRep").innerText = "Por favor, ingresa valores válidos.";
        return;
    }

    let repeticiones = repeticionesTexto.split(",").map(num => parseInt(num.trim()));
    if (repeticiones.some(isNaN)) {
        document.getElementById("resultadoPermRep").innerText = "Las repeticiones deben ser números separados por comas.";
        return;
    }

    let denominador = 1;
    repeticiones.forEach(r => denominador *= factorial(r));

    let resultado = factorial(n) / denominador;
    document.getElementById("resultadoPermRep").innerText = `Resultado: ${resultado}`;
    document.getElementById("interpretacionPermRep").innerText = 
        `Existen ${resultado} formas diferentes de ordenar ${n} elementos con las repeticiones indicadas.`;
}


function calcularProbabilidad() {
    const favorables = parseFloat(document.getElementById("favorable").value);
    const total = parseFloat(document.getElementById("total").value);
    const resultado = document.getElementById("resultadoProb");
    const interpretacion = document.getElementById("interpretacionProb");

    if (isNaN(favorables) || isNaN(total) || total <= 0 || favorables > total) {
        resultado.textContent = "⚠️ Ingresa valores válidos.";
        interpretacion.textContent = "";
        return;
    }

    const prob = favorables / total;
    resultado.textContent = `P(E) = ${prob.toFixed(3)}`;
    interpretacion.textContent = `Esto significa que la probabilidad del evento es del ${(prob * 100).toFixed(2)}%.`;
}



function calcularCondicional() {
    let inter = parseFloat(document.getElementById("interseccion").value);
    let pb = parseFloat(document.getElementById("eventoB").value);

    if (isNaN(inter) || isNaN(pb) || inter < 0 || pb <= 0 || inter > pb) {
        document.getElementById("resultadoCondicional").innerText = 
            "Valores inválidos. Asegúrate que P(A ∩ B) ≤ P(B).";
        return;
    }

    let resultado = inter / pb;
    document.getElementById("resultadoCondicional").innerText =
        `Resultado: P(A | B) = ${resultado.toFixed(3)}`;
}


// ----- TEOREMA DE BAYES -----
function calcularBayes() {
    let pA = parseFloat(document.getElementById("pA").value);
    let pBgA = parseFloat(document.getElementById("pBgivenA").value);
    let pBgNotA = parseFloat(document.getElementById("pBgivenNotA").value);

    if (isNaN(pA) || isNaN(pBgA) || isNaN(pBgNotA) ||
        pA < 0 || pA > 1 || pBgA < 0 || pBgA > 1 || pBgNotA < 0 || pBgNotA > 1) {
        document.getElementById("resultadoBayes").innerText = "Ingresa valores válidos entre 0 y 1.";
        return;
    }

    let pNotA = 1 - pA;

    // Calcular P(B)
    let pB = (pBgA * pA) + (pBgNotA * pNotA);

    // Calcular Bayes
    let resultado = (pBgA * pA) / pB;

    document.getElementById("resultadoBayes").innerText = `P(A|B) = ${resultado.toFixed(4)}`;
    document.getElementById("interpretacionBayes").innerText =
        `La probabilidad actualizada de que ocurra A dado que ocurrió B es del ${(resultado * 100).toFixed(2)}%.`;
}



// =====================================
// FUNCIONES AUXILIARES
// =====================================

// factorial simple
function fact(n) {
    if (n < 0) return 0;
    let res = 1;
    for (let i = 1; i <= n; i++) res *= i;
    return res;
}

// combinatoria nCx
function comb(n, x) {
    if (x < 0 || x > n) return 0;
    return fact(n) / (fact(x) * fact(n - x));
}



// =====================================
// A) DISTRIBUCIÓN BINOMIAL
// =====================================

function calcularBinomial() {
    const n = parseInt(document.getElementById("bin_n").value);
    const p = parseFloat(document.getElementById("bin_p").value);
    const x = parseInt(document.getElementById("bin_x").value);
    const tipo = document.getElementById("bin_tipo").value;

    let resultado = 0;

    function px(k) {
        return comb(n, k) * Math.pow(p, k) * Math.pow(1 - p, n - k);
    }

    if (tipo === "exacto") {
        resultado = px(x);
    }

    if (tipo === "al_menos") {
        for (let k = x; k <= n; k++) {
            resultado += px(k);
        }
    }

    if (tipo === "a_lo_sumo") {
        for (let k = 0; k <= x; k++) {
            resultado += px(k);
        }
    }

    document.getElementById("bin_res").textContent =
        "Resultado: " + resultado.toFixed(6);
}



// =====================================
// B) DISTRIBUCIÓN BINOMIAL NEGATIVA
// =====================================
//
// P(K = k) = C(k-1, r-1) * p^r * (1-p)^(k-r)

function calcularBinomialNegativa() {
    const p = parseFloat(document.getElementById("bn_p").value);
    const r = parseInt(document.getElementById("bn_r").value);
    const k = parseInt(document.getElementById("bn_k").value);
    const tipo = document.getElementById("bn_tipo").value;

    let resultado = 0;

    function prob(k) {
        if (k < r) return 0;
        return comb(k - 1, r - 1) * Math.pow(p, r) * Math.pow(1 - p, k - r);
    }

    if (tipo === "exacto") {
        resultado = prob(k);
    }

    if (tipo === "a_lo_sumo") {
        for (let i = r; i <= k; i++) resultado += prob(i);
    }

    if (tipo === "al_menos") {
        for (let i = k; i <= k + 100; i++) resultado += prob(i);
    }

    document.getElementById("bn_res").textContent =
        "Resultado: " + resultado.toFixed(6);
}



// =====================================
// C) DISTRIBUCIÓN DE POISSON
// =====================================
//
// P(X = k) = e^-λ * λ^k / k!

function calcularPoisson() {
    const lam = parseFloat(document.getElementById("po_lam").value);
    const k = parseInt(document.getElementById("po_k").value);
    const tipo = document.getElementById("po_tipo").value;

    function px(k) {
        return Math.exp(-lam) * Math.pow(lam, k) / fact(k);
    }

    let resultado = 0;

    if (tipo === "exacto") {
        resultado = px(k);
    }

    if (tipo === "menos") {
        for (let i = 0; i < k; i++) resultado += px(i);
    }

    if (tipo === "mas") {
        for (let i = k + 1; i <= k + 100; i++) resultado += px(i);
    }

    document.getElementById("po_res").textContent =
        "Resultado: " + resultado.toFixed(6);
}


// =============================
// EXPONENCIAL
// =============================
function mostrarEntreExponencial() {
    document.getElementById("exp_extra").style.display =
        document.getElementById("exp_tipo").value === "entre" ? "block" : "none";
}

function calcularExponencial() {
    const λ = parseFloat(document.getElementById("exp_lam").value);
    const t = parseFloat(document.getElementById("exp_t").value);
    const tipo = document.getElementById("exp_tipo").value;

    if (isNaN(λ) || λ <= 0) {
        document.getElementById("exp_res").textContent = "λ debe ser mayor que 0.";
        return;
    }

    let res = 0;

    if (tipo === "menor") res = 1 - Math.exp(-λ * t);
    if (tipo === "mayor") res = Math.exp(-λ * t);

    if (tipo === "entre") {
        const a = parseFloat(document.getElementById("exp_a").value);
        const b = parseFloat(document.getElementById("exp_b").value);
        res = Math.exp(-λ * a) - Math.exp(-λ * b);
    }

    document.getElementById("exp_res").textContent = "Resultado: " + res.toFixed(6);
}


// =============================
// NORMAL (usa erf)
// =============================

// CDF de la normal estándar
function normalCDF(z) {
    return 0.5 * (1 + erf(z / Math.sqrt(2)));
}

// función erf (aproximación)
function erf(x) {
    const a1 = 0.254829592,
          a2 = -0.284496736,
          a3 = 1.421413741,
          a4 = -1.453152027,
          a5 = 1.061405429,
          p = 0.3275911;

    const sign = x < 0 ? -1 : 1;
    x = Math.abs(x);
    const t = 1 / (1 + p * x);

    const y = 1 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*Math.exp(-x*x);

    return sign * y;
}

function mostrarEntreNormal() {
    const tipo = document.getElementById("norm_tipo").value;
    document.getElementById("norm_extra").style.display = tipo === "entre" ? "block" : "none";
    document.getElementById("norm_inv").style.display = tipo === "inverso" ? "block" : "none";
}

function calcularNormal() {
    const μ = parseFloat(document.getElementById("norm_mu").value);
    const σ = parseFloat(document.getElementById("norm_sigma").value);
    const x = parseFloat(document.getElementById("norm_x").value);
    const tipo = document.getElementById("norm_tipo").value;

    if (σ <= 0) {
        document.getElementById("norm_res").textContent = "σ debe ser mayor que 0.";
        return;
    }

    let z, res;

    if (tipo === "menor") {
        z = (x - μ) / σ;
        res = normalCDF(z);
    }

    if (tipo === "mayor") {
        z = (x - μ) / σ;
        res = 1 - normalCDF(z);
    }

    if (tipo === "entre") {
        const a = parseFloat(document.getElementById("norm_a").value);
        const b = parseFloat(document.getElementById("norm_b").value);
        const za = (a - μ) / σ;
        const zb = (b - μ) / σ;
        res = normalCDF(zb) - normalCDF(za);
    }

    if (tipo === "inverso") {
        const p = parseFloat(document.getElementById("norm_prob").value);
        z = Math.sqrt(2) * inverseErf(2 * p - 1);
        res = μ + z * σ;
    }

    document.getElementById("norm_res").textContent = "Resultado: " + res.toFixed(6);
}

// inversa de erf (para Z inverso)
function inverseErf(x) {
    const a = 0.147;
    const ln = Math.log((1 - x)*(1 + x));
    return Math.sign(x) * Math.sqrt( Math.sqrt((2/(Math.PI*a) + ln/2)**2 - ln/a) - (2/(Math.PI*a) + ln/2) );
}
