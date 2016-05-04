/* Returns a function that takes in xmin, xmax range and creates pdf*/
function normal_pdf(mean, variance) {
	var coeff = 1 / (Math.sqrt(variance) * Math.sqrt(2*Math.PI))
	var denom = 2*variance

	return function pdf_func(xmin, xmax, samples) {
		result = [] // each  item in form of {x: number, y: number}
		var increments = (xmax-xmin) / samples

		for (var i = xmin; i < xmax; i+=increments) {
			var pdf = coeff * Math.exp( -Math.pow(i-mean, 2) / denom )
			result.push({x: i, y: pdf})
		}
		return result
	}
}

/* given data, samples # of points */
function sample(dist, size) {
	var samples = []
	for (var i = 0; i < size; i++) {
		var index = Math.floor(Math.random() * dist.length)
		samples.push(dist[index])
	}
	return samples
}

/* copies the data points and accumulates */
function accumulate(data) {	
	var copied = [data[0]]

	for (var i = 1; i < data.length; i++) {
		var current = data[i]
		var prev = copied[i-1]
		copied.push({x: i, y: current.y + prev.y})
	}

	return copied
}

/* parts of the equation needed to model gumbel distribution */
function risk(u,p,tick,alpha,tau) {
	return -alpha*(Math.exp(tau*(p[tick-1].y - u[tick-1].y)))
}

function momentum(u,p,tick,beta) {
	return beta*(p[tick-1].y - p[tick-2].y)
}

function underlying(u,tick,gamma) {
	return gamma*(u[tick].y-u[tick-1].y)
}

function gumbel_pdf(u, p, tick, alpha, tau, beta, gamma, delta) {
	var coeff = 1 / delta
	var r = risk(u, p, tick, alpha, tau)
	var m = momentum(u, p, tick, beta)
	var u = underlying(u, tick, gamma)
	var q = r+m+u

	return function pdf_func(xmin, xmax, samples) {
		result = [] // each  item in form of {x: number, y: number}
		var increments = (xmax-xmin) / samples

		for (var i = xmin; i < xmax; i+=increments) {
			var z = (i - q) / delta
			var pdf = coeff * Math.exp(-(z + Math.exp(-z)))
			result.push({x: i, y: pdf})
		}
		return result
	}
}

// returns a function that returns a list of cdf points for each x value
function gumbel_cdf(u, p, tick, alpha, tau, beta, gamma, delta) {
	var r = risk(u, p, tick, alpha, tau)
	var m = momentum(u, p, tick, beta)
	var u = underlying(u, tick, gamma)
	var q = r+m+u

	return function cdf_func(xmin, xmax, samples) {
		result = [] // each  item in form of {x: number, y: number}
		var increments = (xmax-xmin) / samples

		for (var i = xmin; i < xmax; i+=increments) {
			var z = (i - q) / delta
			var cdf = Math.exp(-Math.exp(-z))
			var rounded = Math.round(cdf * 10000) / 10000
			result.push({x: i, y: rounded})
		}
		return result
	}
}

/* given two pdf distributions, where the CDF of seller equals 1-CDF of buyer*/
function find_intersect(buyer_params, seller_params, u, p, tick, xmin, xmax, samples) {
	var b_alpha = buyer_params[0] 
	var b_tau = buyer_params[1]
	var b_beta = buyer_params[2]
	var b_gamma = buyer_params[3]
	var b_delta = buyer_params[4]

	var s_alpha = seller_params[0] 
	var s_tau = seller_params[1]
	var s_beta = seller_params[2]
	var s_gamma = seller_params[3]
	var s_delta = seller_params[4]

	var mindist = 5000;
	var minindex = -1;
	var ratio = -1;

	var b_ratio_last = 0
	var s_ratio_last = 0

	var increments = (xmax-xmin) / samples
	for (var i = xmin; i < xmax; i+=increments) {

		var b_r = risk(u, p, tick, b_alpha, b_tau)
		var b_m = momentum(u, p, tick, b_beta)
		var b_u = underlying(u, tick, b_gamma)
		var b_q = b_r+b_m+b_u

		var b_z = (i - b_q) / b_delta
		var b_cdf = Math.exp(-Math.exp(-b_z))
		var b_rounded = Math.round(b_cdf * 10000) / 10000
		b_ratio_last = 1-b_rounded

		var s_r = risk(u, p, tick, s_alpha, s_tau)
		var s_m = momentum(u, p, tick, s_beta)
		var s_u = underlying(u, tick, s_gamma)
		var s_q = s_r+s_m+s_u

		var s_z = (i - s_q) / s_delta
		var s_cdf = Math.exp(-Math.exp(-s_z))
		var s_rounded = Math.round(s_cdf * 10000) / 10000
		s_ratio_last = s_rounded

		var diff = Math.abs((1-b_rounded) - s_rounded)
		if (diff < mindist) {
			mindist = diff
			minindex = i
			ratio = s_rounded
		}
	}

	return [minindex, ratio]
}

