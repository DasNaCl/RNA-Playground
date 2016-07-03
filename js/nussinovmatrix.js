/**
 * @file Main file containing backend algorithms for RNA-algorithms-JS project.
 Main Items contains in theis file are:
 -Matrix class
 -Nussinov algorithms
 -Traceback algorithms
 * @authors "Martin Mann", "Syed Mohsin Ali".
 */
"use strict";

/**
 * Utility class that covers RNA specific functions.
 *
 */
var RnaUtil = {

    /**
     * checks base pair complementary of two nucleotides
     * @param {string} nt1 first nucleotide
     * @param {string} nt2 second nucleotide
     * @returns {boolean} true if complementary; false otherwise
     */
    areComplementary: function (nt1, nt2) {
        //console.log("areComp:", nt1, nt2);
        var complementary =
            (nt1 === "A" && nt2 === "U") || (nt1 === "U" && nt2 === "A") ||
            (nt1 === "G" && nt2 === "C") || (nt1 === "C" && nt2 === "G") ||
            (nt1 === "G" && nt2 === "U") || (nt1 === "U" && nt2 === "G");
        return complementary;
    },

    /**
     * checks whether or not the sequence is a valid RNA sequence
     * @param sequence
     */
    isRnaSequence: function (sequence) {
        var isValid =
            // check if sequence given
            (sequence !== null)
                // check RNA alphabet
            && sequence.match("^[ACGU]+$");
        return isValid;
    }


};


/**
 * Ancestor information for a certain traceback
 */
var NussinovCellTrace = {

// list of parent cells
    parents: null,

// list of base pairs added
    bps: null,

    /**
     * initializes the object
     * @param {object} parents the parents to set
     * @param {object} bps the base pairs to set
     * @returns {NussinovCellTrace} this for call chaining
     */
    init: function (parents, bps) {
        this.parents = null;
        this.bps = null;
        // input check
        if ((parents !== null && bps === null) || (parents === null && bps !== null)) {
            console.log("ERROR : NussinovCellTrace.init : only one value null");
            return this;
        }
        // store sane data
        this.parents = parents;
        this.bps = bps;
        return this;
    }


};

/**
 * data stored within a cell of a nussinov matrix
 */
var NussinovCell = {

// row
    i: -1,

// column
    j: -1,

// value
    value: null,

// traces for the current value
    traces: null,

    /**
     * inits a cell with the given data and sets traces to an empty list
     * @param i the row of the cell
     * @param j the column of the cell
     * @param value the value of the cell
     * @return this : cell access for chaining
     */
    init: function (i, j, value) {
        // init data
        this.i = i;
        this.j = j;
        this.value = value;
        this.traces = [];
        // this access for chaining
        return this;
    }

};

/**
 * Encodes a full/partial traceback for a given cell
 */
var Traceback = {

    /** structure in dot-bracket-notation */
    structure: "",
    /** list of cell traces of the form [source, parent1, parent2, ...] */
    traces: [],
    //partialStructure: [],

    /**
     * set structure value
     */
    setStructure: function (struct) {
        this.structure = struct;
    },

    /**
     * set traces
     */
    setTraces: function (args) {
        if (this.traces == [])
            this.traces = args;
        else this.traces.push(args);
    },

    /**
     * initialize structure with dot string of sequence length
     * @param (int) len length of string to be created.
     */
    init: function (len) {
        if (len <= 0) return "";
        for (var i = 0; i < len; i++) {
            if (this.structure == undefined)
                this.structure = ".";
            else
                this.structure = this.structure + "."
        }
        this.traces = [];
        return this;
    },

    /**
     * Add Brackets to structure for given basepair
     * @param {array} basepair
     */
    addBrackets: function (basepair) {
        if (basepair[1] > this.structure.length || basepair == null) {
            console.log("String too short or empty bp");
            return;
        }

        this.structure = this.structure.substr(0, basepair[0] - 1) + '(' + this.structure.substr(basepair[0], basepair[1] - basepair[0] - 1) + ')' + this.structure.substr(basepair[1]);
        return this.structure;
    }

    //setPartialStructure: function (struct) {
    //
    //    if (this.partialStructure == []) {
    //
    //        this.partialStructure = [struct];
    //    }
    //
    //    else this.partialStructure.push(struct);
    //},

};

/**
 * Nussinov matrix object
 */
var NussinovMatrix = {

        /**
         * Access to the sequence for this matrix
         */
        sequence: null,

        /**
         * Access name of recursion used
         */
        name: null,

        /**
         * Minimal loop length within computation
         */
        minLoopLength: 0,

        /**
         * cells of the matrix
         */
        cells: [],

        /**
         * The latex representation of the formula computing the matrix.
         */
        latex_representation: "$$",

        /**
         * The table description
         */
        descripion: "default description",

        /**
         * Tracebacks allowance
         */
        allowTraceBack: true,


        /**
         * initialize a matrix of dim = n+1 with indices 0..n, where n is the
         * length of the provided sequence
         * @param {string} sequence the RNA sequence (not null or empty)
         * @returns {NussinovMatrix} this
         */
        init: function
            (sequence, name) { //initialize matrix

            // reset data
            this.sequence = null;
            this.name =null;
            this.cells = [];

            // check input
            if (sequence == null || sequence === "" || name == null) {
                console.log("Matrix init failed for sequence (", sequence, ")");
                return this;
            }

            // store sequence
            this.sequence = sequence;
            this.name = name;

            // create matrix cells
            var n = this.sequence.length;
            for (var i = 0; i <= n; i++) {
                this.cells[i] = [];
                for (var j = 0; j <= n; j++) {
                    // create new cell and initialize
                    if(this.name === "structuresCount") {
                        this.cells[i][j] = Object.create(NussinovCell).init(i, j, null);
                    }
                    else if(this.name === "McKaskill" || this.name === "McKaskill Base"){
                        this.cells[i][j] = Object.create(NussinovCell).init(i, j, null);
                    }
                    else {
                        this.cells[i][j] = Object.create(NussinovCell).init(i, j, 0);
                    }
                }
                ;
            }
            ;

            return this;
        }
        ,

        /**
         * access whole object at cell location (i,j) in matrix
         * @param {int} i row #.
         * @param {int} j column #.
         * @returns {NussinovCell} cell object or null if not available
         */
        getCell: function (i, j) {
            // check border cases
            if (i < 0 || j < 0 || i >= this.getDim() || j >= this.getDim()) {
                return null;
            }
            // check invalid cell access
            if (i > j + 1) {
                return null;
            }
            return this.cells[i][j];
        }
        ,

        /**
         * Gives the dimension of the quadratic matrix
         * @returns {int} the dimension (#rows and columns)
         */
        getDim: function () {
            if (this.cells === null) {
                return 0;
            }
            return this.cells.length;
        }
        ,


        /**
         * returns traces of cell (i,j), i.e. ancestor cells and the basepairs they form
         * @param {int} i row #.
         * @param {int} j column #.
         * @returns {object} ancestor's object Eg. {parents:[],bPs:[]} or null if not available
         */
        getTraces: function (i, j) { //get traceback info for each cell
            // access cell at location (i,j) in the matrix
            var cell = this.getCell(i, j);
            // check if valid cell
            if (cell === null) {
                return null;
            }
            return cell.traces;
        }
        ,

        computeValue: function(i, j) {
            // Computes the M[i, j] by accessing cells using getValue function
            // for memoization.
            return 0;
        },
        /**
         * access value of cell at location (i,j) in matrix
         * @param {int} i row #.
         * @param {int} j column #.
         * @returns {int} cell value or 0 if invalid cell
         */
        getValue: function (i, j) {
            // access cell at location (i,j) in the matrix
            var cell = this.getCell(i, j);
            if (cell === null) {
                return null;
            }
            // check if invalid cell
            if (cell.value === null) {
                cell.value = this.computeValue(i, j);
            }
            // get cell value
            return cell.value;
        }
        ,


        /**
         * Updates the ancestor list of a given cell if the curVal is higher or
         * equal to the current value within the cell.
         * If the value is equal, curAncestor is added to the list.
         * If the value is smaller than curVal, curAncestor will be set to be the
         * only list entry.
         * @param i the row of the cell to update
         * @param j the column of the cell to update
         * @param curAncestor
         */
        updateCell: function (i, j, curAncestor) {
            // get cell to update
            var curCell = this.getCell(i, j);
            // check if something to update
            if (curCell === null) {
                return;
            }

            // init value with number of additional base pairs
            var curVal = curAncestor.bps.length;
            // add scores of ancestor cells
            for (var x = 0; x < curAncestor.parents.length; x++) {
                curVal += this.getValue(curAncestor.parents[x][0], curAncestor.parents[x][1]);
            }
            // check if we have to update
            if (curCell.value <= curVal) {
                // check for new maximal value
                if (curCell.value < curVal) {
                    // reset ancestor list
                    curCell.traces = [];
                    // store new maximum
                    curCell.value = curVal;
                }
                ;
                // store this ancestor
                curCell.traces.push(curAncestor);
            }
            ;

        }
        ,


        /**
         * Fills the matrix according to the recursion.
         *
         * NOTE: this function has to be overwritten by subclasses
         *
         * @param {string} sequence the RNA sequence to compute the matrix for
         * @param {int} minLoopLength the minimal loop length to be used for computation
         *
         * @returns {NussinovMatrix} this for call chaining
         */
        computeMatrix: function (sequence, minLoopLength) {
            console.log("WARNING: computeMatrix() not implemented in NussinovMatrix superclass; overwrite in subclass!");

            // resize and init matrix
            this.init(sequence);

            // set minimal loop length
            this.minLoopLength = minLoopLength;

            return this;
        }
        ,


        /**
         * Access to the recursion in LaTeX encoding that is used in the computeMatrix implementation
         *
         * NOTE: this function has to be overwritten by subclasses
         *
         * @returns {string} latex encoding of the recursion
         */
        getRecursionInLatex: function () {
            console.log("WARNING: getRecursionInLatex() not implemented in NussinovMatrix superclass; overwrite in subclass!");
            return "";
        }
        ,


        /**
         * Returns a description for the implemented recursion
         *
         * NOTE: this function has to be overwritten by subclasses
         *
         * @returns {string} description of the recursion
         */
        getRecursionDescription: function () {
            console.log("WARNING: getRecursionDescription() not implemented in NussinovMatrix superclass; overwrite in subclass!");
            return "";
        }
        ,

        /**
         * creates a string representation of the matrix
         * @returns {string} matrix as string
         */
        toString: function () {
            var str = this.minLoopLength + "    ";
            for (var i = 0; i + 1 < this.getDim(); i++) {
                str += this.sequence[i] + "  ";
            }
            str += "\n";
            for (var i = 1; i < this.getDim(); i++) {
                // print sequence
                if (i === 0) {
                    str += "  ";
                } else {
                    str += this.sequence[i - 1] + " ";
                }
                // print values
                for (var j = 0; j < this.getDim(); j++) {
                    if (j !== 0) {
                        str += ", ";
                    }
                    // get value
                    str += this.getValue(i, j);
                    //    			/* get number of traces */
                    //    			if (this.getCell(i,j) === null || this.getCell(i,j).traces === null) {
                    //    				str += "0";
                    //    			} else {
                    //    				str += this.getCell(i,j).traces.length;
                    //    			}
                }
                ;
                str += "\n";
            }
            ;
            return str;
        }
        ,


        /**
         * Returns ONE optimal traceback if the matrix was already computed.
         *
         * @returns a Traceback object or null if the matrix wasnt computed yet
         */
        getOptimalTraceback: function (i, j, optimalTrace) {
            // check if there was something computed yet
            if (this.sequence === null) {
                return null;
            }
            ;

            // create empty trace to fill
            if (typeof optimalTrace == "undefined") {
                //console.log("optimalTrace object is not defined");
            }
            ;

            if (this.getValue(i, j) === null) {
                return;
            }
            else {
                //console.log(this.getTraces(i,j)[0]);

                for (var p in this.getTraces(i, j)[0].parents) {
                    var parent = this.getTraces(i, j)[0].parents[p];
                    var basepair = this.getTraces(i, j)[0].bps;

                    if (typeof basepair[0] != "undefined") {
                        optimalTrace.addBrackets(basepair[0]);
                    }
                    ;

                    this.getOptimalTraceback(parent[0], parent[1], optimalTrace);

                }
                ;
                optimalTrace.traces.push([[i, j], this.getTraces(i, j)[0].parents]);

            }
            ;

            return optimalTrace;
        }
        ,
        /**
         * countBasepairs(for wuchty)
         */
        countBasepairs: function(bps, sigma){
            var NSprime = bps.length;
            for (var s in sigma) {
                var i = sigma[s][0];
                var j = sigma[s][1];

                NSprime += this.getValue(i, j);
            }
            return NSprime;
        },




        /**
         * WUCHTY
         */
        wuchty: function (delta) {
            var seq_length = this.sequence.length;
            var Nmax = this.getCell(1, seq_length).value;

            var S = {sigma: [[1, seq_length]], P: [], traces: []};
            var R = [S];
            var SOS = [];
            var loop = 0;

            while (R.length != 0) {

                // Pop R
                var pop_R = R.pop();
                var sigma = pop_R.sigma;
                var P = pop_R.P;
                var t_traces = JSON.stringify(pop_R.traces);
                var traces = JSON.parse(t_traces);

                var sigma_remaining = 0;
                for (var s in sigma) {//console.log("var s:", sigma[s]);
                    if ((sigma[s][0]) <= (sigma[s][1] - this.minLoopLength)) sigma_remaining++;
                }

                if (sigma.length == 0 || sigma_remaining == 0) {
                    var temp_sos = {structure: this.conv_str(P, seq_length), traces: traces};
                    SOS.push(temp_sos);
                }

                else {
                    var ij = sigma.pop();

                    if (ij[0] > ij[1]) {
                        pop_R.traces = traces;
                        R.push(pop_R);
                        continue;
                    }

                    var json_ij_traces = JSON.stringify(this.getCell(ij[0], ij[1]).traces);
                    var ij_traces = JSON.parse(json_ij_traces);

                    var S_prime = {};

                    for (var t in ij_traces) {
                        var S_prime = {};
                        var t_trace = JSON.stringify(ij_traces[t]);
                        var trace = JSON.parse(t_trace);
                        var trace_p = JSON.parse(t_trace);

                        // add sigma info in S_prime
                        S_prime.sigma = [];
                        if (sigma.length > 0) {
                            for (var s in sigma)
                                trace.parents.unshift(sigma[s]);
                            S_prime.sigma = trace.parents;
                        }
                        else S_prime.sigma = trace.parents;

                        // add P(bps) info in S_prime
                        S_prime.P = []
                        if (P[0] != undefined) {
                            for (var p in P) {
                                S_prime.P.push(P[p]);
                            }
                        }
                        if (trace.bps[0] != undefined) {
                            S_prime.P.push(trace.bps[0]);
                        }

                        // add traces info in S_prime
                        var temp_trace = [ij];
                        temp_trace.push(trace_p.parents);
                        if (traces.length == 0) {
                            S_prime.traces = [temp_trace];
                        }
                        else {
                            var clone_traces = JSON.stringify(traces);
                            var parse_clone_traces = JSON.parse(clone_traces);
                            parse_clone_traces.unshift(temp_trace);
                            S_prime.traces = parse_clone_traces;
                        }
                        R.push(S_prime);
                    }
                }
                ;

            }
            return SOS;
        },

        conv_str: function (x, length) {
            var str = "";
            for (var l = 0; l < length; l++) {
                str += ".";
            }
            //console.log(str, length);
            for (var i in x) {
                str = str.substr(0, x[i][0] - 1) + "(" + str.substr(x[i][0], str.length - 1);
                str = str.substr(0, x[i][1] - 1) + ")" + str.substr(x[i][1], str.length - 1);
            }
            return str;
        },

    }
    ;

var DPAlgorithm = {
    Description: "Algorithm",

    //Tables: [],

    defaultPars: {},

    UpdateCells: function(i, j){},

    computeMatrices: function(sequence, args){},
};

/****** NussinovMatrix_ambiguous extending NussinovMatrix ************************/

/**
 * Implements an ambiguous recursion
 */
var NussinovMatrix_ambiguous = Object.create(NussinovMatrix);


/**********  OVERWRITING + EXTENSION  ********************/

/**
 * Returns a description for the implemented recursion
 *
 * @returns {string} description of the recursion
 */
NussinovMatrix_ambiguous.getRecursionDescription = function () {
    return "Ambiguous recursion";
};


/**
 * Access to the recursion in LaTeX encoding that is used in the computeMatrix implementation
 *
 * @returns {string} latex encoding of the recursion
 */
NussinovMatrix_ambiguous.getRecursionInLatex = function () {
    return "$D(i,j) = \\max \\begin{cases} D(i+1,j) & S_i \\text{ unpaired} \\\\ D(i,j-1) & S_j \\text{ unpaired} \\\\ D(i+1,j-1)+1 &  S_i,S_j \\text{ compl. base pair and } i+ l< j \\\\ \\max_{i< k< (j-1)} D(i,k)+D(k+1,j) & \\text{decomposition} \\end{cases}$";
};

/**
 * Fills the matrix according to the recursion.
 *
 * @param {string} sequence the RNA sequence to compute the matrix for
 * @param {int} minLoopLength the minimal loop length to be used for computation
 *
 * @returns {NussinovMatrix} this for call chaining
 */
NussinovMatrix_ambiguous.computeMatrix = function (sequence, minLoopLength) {
// resize and initialize matrix
    this.init(sequence, "ambiguous");
// store minimal loop length
    this.minLoopLength = minLoopLength;
    //console.log("computing ambiguos matrix");
// fill matrix by diagonals
// iterate over all substructure spans that can have a base pair
    for (var span = minLoopLength; span < this.getDim(); span++) {
        // iterate over all rows
        for (var i = 1; i < this.getDim() - minLoopLength; i++) {
            // get column for current span
            var j = i + span;

            // i unpaired
            this.updateCell(i, j, Object.create(NussinovCellTrace).init([[i + 1, j]], []));

            // j unpaired
            this.updateCell(i, j, Object.create(NussinovCellTrace).init([[i, j - 1]], []));

            // check (i,j) base pair
            if ((j - i > minLoopLength) && RnaUtil.areComplementary(this.sequence[i - 1], this.sequence[j - 1])) {
                // get value for base pair
                this.updateCell(i, j, Object.create(NussinovCellTrace).init([[i + 1, j - 1]], [[i, j]]));
            }
            ;

            // check decomposition into substructures (minLength==2)
            for (var k = i + 1; k < (j - 1); k++) {
                // get decomposition value
                this.updateCell(i, j, Object.create(NussinovCellTrace).init([[i, k], [k + 1, j]], []));
            }
            ;
        }
        ;
    }
    ;

    return this;
};

/**
 * getSubstructures(for wuchty)
 */
NussinovMatrix_ambiguous.getSubstructures = function(sigma, P, traces, delta, maxLengthR ) {
    var Nmax = this.getCell(1, this.sequence.length).value;
    var R = [];
    var ij = sigma.pop();
    //console.log(ij);

    // check for sane interval
    // if i>j dont continue
    if (ij[0] >= ij[1]) {
        //console.log("ij[0] > ij[1]", ij[0], ij[1]);
        var S_prime = {};

        var sigma_prime = JSON.stringify(sigma);
        sigma_prime = JSON.parse(sigma_prime);

        var tmp_P = JSON.stringify(P);
        tmp_P = JSON.parse(tmp_P);

        var tmp_traces = JSON.stringify(traces);
        tmp_traces = JSON.parse(tmp_traces);

        S_prime.sigma = sigma_prime;
        S_prime.P = tmp_P;
        S_prime.traces = tmp_traces;

        R.push(S_prime);
        //console.log("returning R:", JSON.stringify(R));
        return R;
    }

    // if (i,j) == (i+1,j-1) + bp(ij)
    {
        if (ij[1] - ij[0] > this.minLoopLength) {
        	//console.log(this.sequence);
        	//console.log(this.sequence[ij[0] - 1], this.sequence[ij[1] - 1]);
            if (RnaUtil.areComplementary(this.sequence[ij[0] - 1], this.sequence[ij[1] - 1])) {
                var sigma_prime = JSON.stringify(sigma);
                sigma_prime = JSON.parse(sigma_prime);
                sigma_prime.push([ij[0] + 1, ij[1] - 1]);

                var tmp_P = JSON.stringify(P);
                tmp_P = JSON.parse(tmp_P);
                tmp_P.push([ij[0], ij[1]]);

                var tmp_traces = JSON.stringify(traces);
                tmp_traces = JSON.parse(tmp_traces);

                var NSprime = this.countBasepairs(tmp_P, sigma_prime);

                if (NSprime >= Nmax - delta) {
                    var S_prime = {};
                    S_prime.sigma = sigma_prime;
                    S_prime.P = tmp_P;
                    tmp_traces.unshift([ij, [[ij[0] + 1, ij[1] - 1]]]);
                    S_prime.traces = tmp_traces;
                    //console.log("i+1,j-1:", JSON.stringify(S_prime));
                    // push to the front to keep base pair most prominent to refine
                    R.unshift(S_prime);
                }
            }
        }

        // check if enough structures found so far
        if (R.length >= maxLengthR) {
        	//console.log("returning R:", JSON.stringify(R));
        	return R;
        }
    }

    // if (i,j) == (i+1,j)
    {
        var sigma_prime = JSON.stringify(sigma);
        sigma_prime = JSON.parse(sigma_prime);
        sigma_prime.unshift([ij[0] + 1, ij[1]]);

        var tmp_P = JSON.stringify(P);
        tmp_P = JSON.parse(tmp_P);

        var tmp_traces = JSON.stringify(traces);
        tmp_traces = JSON.parse(tmp_traces);

        var NSprime = this.countBasepairs(P, sigma_prime);

        if (NSprime >= Nmax - delta) {
            var S_prime = {};
            S_prime.sigma = sigma_prime;
            S_prime.P = tmp_P;
            tmp_traces.unshift([ij, [[ij[0] + 1, ij[1]]]]);
            S_prime.traces = tmp_traces;
            //console.log("i+1,j:", JSON.stringify(S_prime));
            // push to the front to keep base pair most prominent to refine
            R.unshift(S_prime);
        }

        // check if enough structures found so far
        if (R.length >= maxLengthR) {
        	//console.log("returning R:", JSON.stringify(R));
        	return R;
        }
    }

    // if (i,j) == (i,j-1)
    {
	    var sigma_prime = JSON.stringify(sigma);
	    sigma_prime = JSON.parse(sigma_prime);
	    sigma_prime.unshift([ij[0], ij[1] - 1]);

	    var tmp_P = JSON.stringify(P);
	    tmp_P = JSON.parse(tmp_P);

	    var tmp_traces = JSON.stringify(traces);
	    tmp_traces = JSON.parse(tmp_traces);

	    var NSprime = this.countBasepairs(P, sigma_prime);

	    if (NSprime >= Nmax - delta) {
	        var S_prime = {};
	        S_prime.sigma = sigma_prime;
	        S_prime.P = tmp_P;
	        tmp_traces.unshift([ij, [[ij[0], ij[1] - 1]]]);
	        S_prime.traces = tmp_traces;
	        //console.log("i,j-1:", JSON.stringify(S_prime));
	        // push to the front to keep base pair most prominent to refine
	        R.unshift(S_prime);
	    }

	    // check if enough structures found so far
	    if (R.length >= maxLengthR) {
        	//console.log("returning R:", JSON.stringify(R));
	    	return R;
	    }
    }

    // if (i,j) == (i,l) + (l+1, j)
    for (var l = ij[0] + 1; l < ij[1] - 1; l++) {

            var sigma_prime = JSON.stringify(sigma);
            sigma_prime = JSON.parse(sigma_prime);
            sigma_prime.push([ij[0], l]);
            sigma_prime.push([l + 1, ij[1]]);

            var tmp_P = JSON.stringify(P);
            tmp_P = JSON.parse(tmp_P);

            var tmp_traces = JSON.stringify(traces);
            tmp_traces = JSON.parse(tmp_traces);

            var NSprime = this.countBasepairs(tmp_P, sigma_prime);

            if (NSprime >= Nmax - delta) {

                var S_prime = {};
                S_prime.sigma = sigma_prime;
                S_prime.P = tmp_P;
                tmp_traces.unshift([ij, [[ij[0], l], [l + 1, ij[1]]]]);
                S_prime.traces = tmp_traces;
                //console.log("ilj:", JSON.stringify(S_prime));
                // push to the front to keep base pair most prominent to refine
                R.unshift(S_prime);
            }

            // check if enough structures found so far
            if (R.length >= maxLengthR) {
            	//console.log("returning R:", JSON.stringify(R));
            	return R;
            }

        }

    //console.log("returning R:", JSON.stringify(R));
    return R;
}
;
/****** NussinovMatrix_unique extending NussinovMatrix ************************/

/**
 * Implements a non-ambiguous recursion
 */
var NussinovMatrix_unique = Object.create(NussinovMatrix);


/**********  OVERWRITING + EXTENSION  ********************/

/**
 * Returns a description for the implemented recursion
 *
 * @returns {string} description of the recursion
 */
NussinovMatrix_unique.getRecursionDescription = function () {
    return "Recursion by Nussinov et al. (1978) with unique decomposition";
};

/**
 * Access to the recursion in LaTeX encoding that is used in the computeMatrix implementation
 *
 * @returns {string} latex encoding of the recursion
 */
NussinovMatrix_unique.getRecursionInLatex = function () {
    return "$D(i,j) = \\max \\begin{cases} D(i,j-1) & S_j \\text{ unpaired} \\\\ \\max_{i\\leq k< (j-l)} D(i,k-1)+D(k+1,j-1)+1 & S_k,S_j \\text{ compl. base pair} \\end{cases}$";
};

/**
 * Fills the matrix according to the recursion.
 *
 * @param {string} sequence the RNA sequence to compute the matrix for
 * @param {int} minLoopLength the minimal loop length to be used for computation
 *
 * @returns {NussinovMatrix} this for call chaining
 */
NussinovMatrix_unique.computeMatrix = function (sequence, minLoopLength) {

// resize and initialize matrix
    this.init(sequence, "unique");
// store minimal loop length
    this.minLoopLength = minLoopLength;

// fill matrix by diagonals
// iterate over all substructure spans that can have a base pair
    for (var span = minLoopLength; span < this.getDim(); span++) {
        // iterate over all rows
        for (var i = 0; i < this.getDim() - minLoopLength; i++) {
            // get column for current span
            var j = i + span;

            // j unpaired
            this.updateCell(i, j, Object.create(NussinovCellTrace).init([[i, j - 1]], []));

            // check base pair based decomposition : (k,j) base pair
            for (var k = i; k + minLoopLength < j; k++) {
                // check if sequence positions are compatible
                if (RnaUtil.areComplementary(this.sequence[k - 1], this.sequence[j - 1])) {
                    this.updateCell(i, j, Object.create(NussinovCellTrace).init([[i, k - 1], [k + 1, j - 1]], [[k, j]]));
                }
                ;
            }
            ;
        }
        ;
    }
    ;

    return this;
};

/**
 * getSubstructures(for wuchty)
 */
NussinovMatrix_unique.getSubstructures = function(sigma, P, traces, delta, maxLengthR) {
    var Nmax = this.getCell(1, this.sequence.length).value;
    var R = [];
    var ij = sigma.pop();

    // if i>j dont countinue
    if (ij[0] > ij[1]) {
        //console.log("ij[0] > ij[1]", ij[0], ij[1]);
        var S_prime = {};
        S_prime.sigma = sigma;
        S_prime.P = P;
        S_prime.traces = traces;
        R.push(S_prime);
        //console.log("returning R:", JSON.stringify(R));
        return R;
    }

    // if (i,j) == (i,l-1) + (l+1, j-1) + 1
    for (var l = ij[0]; l <= ij[1] - 1; l++) {
        if (ij[1] - l > this.minLoopLength) {
            if (RnaUtil.areComplementary(this.sequence[l - 1], this.sequence[ij[1] - 1])) {
                var sigma_prime = JSON.stringify(sigma);
                sigma_prime = JSON.parse(sigma_prime);
                sigma_prime.push([ij[0], l - 1]);
                sigma_prime.push([l + 1, ij[1] - 1]);

                var tmp_P = JSON.stringify(P);
                tmp_P = JSON.parse(tmp_P);
                tmp_P.push([l, ij[1]]);

                var tmp_traces = JSON.stringify(traces);
                tmp_traces = JSON.parse(tmp_traces);

                var NSprime = this.countBasepairs(tmp_P, sigma_prime);

                if (NSprime >= Nmax - delta) {

                    var S_prime = {};
                    S_prime.sigma = sigma_prime;
                    S_prime.P = tmp_P;
                    tmp_traces.unshift([ij, [[ij[0], l - 1], [l + 1, ij[1] - 1]]]);
                    S_prime.traces = tmp_traces;
                    //console.log("ilj:", JSON.stringify(S_prime));
                    // push to the back to keep base pair most prominent to refine
                    R.push(S_prime);
                }
                
                // check if enough structures found so far
                if (R.length >= maxLengthR) {
                    //console.log("returning R:", JSON.stringify(R));
                	return R;
                }

            }
        }
    }
    
    // if (i,j) == (i,j-1)
    {
	    var sigma_prime = JSON.stringify(sigma);
	    sigma_prime = JSON.parse(sigma_prime);
	    sigma_prime.unshift([ij[0], ij[1] - 1]);
	
	    var tmp_P = JSON.stringify(P);
	    tmp_P = JSON.parse(tmp_P);
	
	    var tmp_traces = JSON.stringify(traces);
	    tmp_traces = JSON.parse(tmp_traces);
	
	    var NSprime = this.countBasepairs(P, sigma_prime);
	
	    if (NSprime >= Nmax - delta) {
	        var S_prime = {};
	        S_prime.sigma = sigma_prime;
	        S_prime.P = tmp_P;
	        tmp_traces.unshift([ij, [[ij[0], ij[1] - 1]]]);
	        S_prime.traces = tmp_traces;
	        //console.log("ij-1:", JSON.stringify(S_prime));
            // push to the front to keep base pair most prominent to refine
	        R.unshift(S_prime);
	    }
	    
	    // check if enough structures found so far
	    if (R.length >= maxLengthR) {
	        //console.log("returning R:", JSON.stringify(R));
	    	return R;
	    }
    }

    //console.log("returning R:", JSON.stringify(R));
    return R;
}
;


/*****************************************************************************/

/**
 * WUCHTY(2nd version) enumerating up to 10 structures
 */
function wuchty_2nd(xmat, delta, formula) {
	// get maximal number of structures to report
	var maxSOS = 10;
	// call subroutine
	return wuchty_2nd_limited(xmat, delta, formula, maxSOS);
}

/**
 * WUCHTY(2nd version)
 */
function wuchty_2nd_limited(xmat, delta, formula, maxSOS) {
    if(xmat.sequence == undefined)return;
    var seq_length = xmat.sequence.length;
    var Nmax = xmat.getCell(1, seq_length).value;
    //console.log("Wuchty beginning\nSequence:", xmat.sequence, "\nNmax:", Nmax, "\nDelta:", delta, "\n");

    var S = {sigma: [[1, seq_length]], P: [], traces: []};
    var R = [S];
    var SOS = [];
    var loop = 0;

    while (R.length != 0) {
        //console.log("\nloop:", ++loop);
        //console.log("R length:" + R.length);

        // Pop R
        var pop_R = R.pop();
        var sigma = pop_R.sigma;
        var P = pop_R.P;
        var t_traces = JSON.stringify(pop_R.traces);
        var traces = JSON.parse(t_traces);
        //console.log("poped R:", JSON.stringify(pop_R));
        //console.log("R_remaining:", JSON.stringify(R));

        var sigma_remaining = 0;
        for (var s in sigma) {//console.log("var s:", sigma[s]);
            if ((sigma[s][0]) <= (sigma[s][1] - xmat.minLoopLength)) sigma_remaining++;
        }

        if (sigma.length == 0 || sigma_remaining == 0) {
            //console.log("no more structs in poped R");
            var temp_sos = {structure: xmat.conv_str(P, seq_length), traces: traces};
            SOS.push(temp_sos);

            //console.log("pushed SOS:", JSON.stringify(SOS));
        }

        else {
            //console.log(formula);
            // compute maximal number of structures still to compute
            var maxLengthR = maxSOS-SOS.length;
            if (maxLengthR < 0 ) maxLengthR = 0;
            var R_prime = formula.getSubstructures(sigma, P, traces, delta, maxLengthR);
            for(var r in R_prime){
                R.push(R_prime[r]);
            }
        }
        // check if enough structures found so far
        if(SOS.length>=maxSOS)
        	break;
        //console.log("R:", JSON.stringify(R));

    }
    //console.log("SOS:", JSON.stringify(SOS));
    //console.log("\nwuchty end");
    //console.log(SOS);
    return SOS;
}
;

/*
 var DPAlgorithm = {
 Description: "Algorithm",

 Tables: [],

 defaultPars: {},

 UpdateCells: function(i, j){},

 computeMatrices: function(sequence, args...){},
 };

 */

/** Nussinov Structures Count*/

var NussinovDPAlgorithm_structuresCount = Object.create(DPAlgorithm);

NussinovDPAlgorithm_structuresCount.Description = "Nussinov counting";

NussinovDPAlgorithm_structuresCount.Tables = new Array();
NussinovDPAlgorithm_structuresCount.Tables.push(Object.create(NussinovMatrix));

NussinovDPAlgorithm_structuresCount.Tables[0].latex_representation = "$C_{i,j} = C_{i,j-1} + \\sum_{i\\leq k <(j-l) \\atop S_k,S_j \\text{ pair}} C_{i,k-1} * C_{k+1,j-1} * 1 $";


NussinovDPAlgorithm_structuresCount.Tables[0].computeValue = function(i, j) {
    if (i > j + 1 || i < 0 || j < 0 || i >= this.getDim() || j >= this.getDim()) {
        return 0;
    }
    if (i >= j) {
        return 1;
    }
    var res = 0;

    console.log(i, j);
    // unpaired
    res += this.getValue(i, j - 1);
    for (var k = i; k < j - this.minLoopLength; ++k) {
        if (RnaUtil.areComplementary(this.sequence[k - 1], this.sequence[j - 1])) {
            res += this.getValue(i, k - 1) * this.getValue(k + 1, j - 1) * 1;
        }
    }

    return res;
};

NussinovDPAlgorithm_structuresCount.computeMatrix = function (sequence, minLoopLength) {
// resize and initialize matrix

    this.Tables[0].init(sequence, "structuresCount");
    // store minimal loop length
    this.Tables[0].minLoopLength = minLoopLength;

    for (var i = 0; i < this.Tables[0].getDim(); i++) {
        for (var j = 0; j < this.Tables[0].getDim(); ++j) {
            // get column for current span
            this.Tables[0].getValue(i, j);
        }
        ;
    }
    ;

    return this.Tables;
};

/*
NussinovMatrix_structuresCount.updateCell = function(i, j, curAncestor){
    var val = 1;
    // get cell to update
    var curCell = this.getCell(i, j);
    // check if something to update
    if (curCell === null) {
        return;
    }

    // add scores of ancestor cells
    for (var x = 0; x < curAncestor.parents.length; x++) {
        if(curAncestor.parents[x][0] > 0 && curAncestor.parents[x][1] > 0
            && curAncestor.parents[x][0] < curAncestor.parents[x][1]) {
            //var q = this.getValue(curAncestor.parents[x][0], curAncestor.parents[x][1]) * Math.exp(-1 * (-1));
            var q = this.getValue(curAncestor.parents[x][0], curAncestor.parents[x][1]);
            val *= q;
        }
    }
    // check for case (i,j-1)
    if(curAncestor.bps.length === 0){
        curCell.value = val;
    }
    else {
        curCell.value += val;
    }
    // curCell.traces.push(curAncestor);
}
;
*/

/****** NussinovMatrix_structuresCount extending NussinovMatrix ************************/

/**
 * Implements a non-ambiguous recursion
 */
var NussinovMatrix_structuresCount = Object.create(NussinovMatrix);


/**********  OVERWRITING + EXTENSION  ********************/

/**
 * Returns a description for the implemented recursion
 *
 * @returns {string} description of the recursion
 */
NussinovMatrix_structuresCount.getRecursionDescription = function () {
    return "Recursion to count number of total structures.";
};

/**
 * Access to the recursion in LaTeX encoding that is used in the computeMatrix implementation
 *
 * @returns {string} latex encoding of the recursion
 */
NussinovMatrix_structuresCount.getRecursionInLatex = function () {
    return "$C_{i,j} = C_{i,j-1} + \\sum_{i\\leq k <(j-l) \\atop S_k,S_j \\text{ pair}} C_{i,k-1} * C_{k+1,j-1} * 1 $";
};

/**
 * Fills the matrix according to the recursion.
 *
 * @param {string} sequence the RNA sequence to compute the matrix for
 * @param {int} minLoopLength the minimal loop length to be used for computation
 *
 * @returns {NussinovMatrix_structuresCount} this for call chaining
 */
NussinovMatrix_structuresCount.computeMatrix = function (sequence, minLoopLength) {

// resize and initialize matrix
    this.init(sequence, "structuresCount");
    //console.log(sequence, this.sequence.length);
// store minimal loop length
    this.minLoopLength = minLoopLength;

// fill matrix by diagonals
// iterate over all substructure spans that can have a base pair
    for (var span = minLoopLength; span < this.getDim(); span++) {
        // iterate over all rows
        for (var i = 0; i < this.getDim() - minLoopLength; i++) {
            // get column for current span
            var j = i + span;

            // j unpaired
            this.updateCell(i, j, Object.create(NussinovCellTrace).init([[i, j - 1]], []));

            // check base pair based decomposition : (k,j) base pair
            for (var k = i; k + minLoopLength < j; k++) {
                // check if sequence positions are compatible
                if (RnaUtil.areComplementary(this.sequence[k - 1], this.sequence[j - 1])) {
                    this.updateCell(i, j, Object.create(NussinovCellTrace).init([[i, k - 1], [k + 1, j - 1]], [[k, j]]));
                }
                ;
            }
            ;
        }
        ;
    }
    ;

    return this;
};
/**
NussinovMatrix_structuresCount.updateCell = function(i, j, curAncestor){
    var val = 1;
    // get cell to update
    var curCell = this.getCell(i, j);
    // check if something to update
    if (curCell === null) {
        return;
    }

    // add scores of ancestor cells
    for (var x = 0; x < curAncestor.parents.length; x++) {
        if(curAncestor.parents[x][0] > 0 && curAncestor.parents[x][1] > 0
            && curAncestor.parents[x][0] < curAncestor.parents[x][1]) {
            //var q = this.getValue(curAncestor.parents[x][0], curAncestor.parents[x][1]) * Math.exp(-1 * (-1));
            var q = this.getValue(curAncestor.parents[x][0], curAncestor.parents[x][1]);
            val *= q;
        }
    }
    // check for case (i,j-1)
    if(curAncestor.bps.length === 0){
        curCell.value = val;
    }
    else {
        curCell.value += val;
    }
    // curCell.traces.push(curAncestor);
}
;
 */


/****** McKaskill_simple extending NussinovMatrix ************************/

var NussinovDPAlgorithm_McKaskill = Object.create(DPAlgorithm);

NussinovDPAlgorithm_McKaskill.Description = "Mckaskill";

NussinovDPAlgorithm_McKaskill.Tables = new Array();
NussinovDPAlgorithm_McKaskill.Tables.push(Object.create(NussinovMatrix));
NussinovDPAlgorithm_McKaskill.Tables.push(Object.create(NussinovMatrix));

NussinovDPAlgorithm_McKaskill.Tables[0].latex_representation = "$$Q_{i,j} = Q_{i,j-1} + \\sum_{i\\leq k <(j-l)} Q_{i,k-1} * Q^{b}_{k,j} $$";
NussinovDPAlgorithm_McKaskill.Tables[1].latex_representation = "$$Q_{i,j}^{b} = \\begin{cases} Q_{i + 1, j - 1} * \\exp(-E(bp)/RT) & \\text{ if }i,j \\text{ can form base pair} \\\\ 0 & \\text{ otherwise}\\end{cases}$$";


NussinovDPAlgorithm_McKaskill.Tables[0].computeValue = function(i, j) {
    if (i > j + 1 || i < 0 || j < 0 || i >= this.getDim() || j >= this.getDim()) {
        return 0;
    }
    if (i >= j) {
        return 1;
    }
    var res = 0;

    console.log(i, j);
    // unpaired
    res += this.getValue(i, j - 1);
    for (var k = i; k < j - this.minLoopLength; ++k) {
        if (RnaUtil.areComplementary(this.sequence[k - 1], this.sequence[j - 1])) {
            res += this.getValue(i, k - 1) * NussinovDPAlgorithm_McKaskill.Tables[1].getValue(k, j);
        }
    }

    return res;
};

NussinovDPAlgorithm_McKaskill.Tables[1].computeValue = function(i, j) {
    if (i > j + 1 || i < 0 || j < 0 || i >= this.getDim() || j >= this.getDim()) {
        return 0;
    }
    if (i >= j) {
        return 1;
    }
    console.log(i, j);
    // unpaired
    return  NussinovDPAlgorithm_McKaskill.Tables[1].getValue(i + 1, j - 1) * Math.exp(1);
};

NussinovDPAlgorithm_McKaskill.computeMatrix = function (sequence, minLoopLength) {

    this.Tables[0].init(sequence, "McKaskill");
    this.Tables[1].init(sequence, "McKaskill Base");
    // store minimal loop length
    this.Tables[0].minLoopLength = minLoopLength;
    this.Tables[1].minLoopLength = minLoopLength;

    for (var i = 0; i < this.Tables[0].getDim(); i++) {
        for (var j = 0; j < this.Tables[0].getDim(); ++j) {
            // get column for current span
            this.Tables[0].getValue(i, j);
        }
        ;
    }
    ;

    for (var i = 0; i < this.Tables[1].getDim(); i++) {
        for (var j = 0; j < this.Tables[1].getDim(); ++j) {
            // get column for current span
            this.Tables[1].getValue(i, j);
        }
        ;
    }
    ;

    return this.Tables;
};
/**
 * Implements a non-ambiguous recursion
 */
var McKaskill_simple = Object.create(NussinovMatrix);

/**
 * Returns a description for the implemented recursion
 *
 * @returns {string} description of the recursion
 */
McKaskill_simple.getRecursionDescription = function () {
    return "Recursion to count the energy of all the structures.";
};

/**
 * Access to the recursion in LaTeX encoding that is used in the computeMatrix implementation
 *
 * @returns {string} latex encoding of the recursion
 */
McKaskill_simple.getRecursionInLatex = function () {
    return "$$Q_{i,j} = Q_{i,j-1} + \\sum_{i\\leq k <(j-l)} Q_{i,k-1} * Q^{b}_{k,j} $$";
};


/**********  OVERWRITING + EXTENSION  ********************/

McKaskill_simple.updateCell = function(i, j, curAncestor){
    var val = 1;
    // get cell to update
    var curCell = this.getCell(i, j);
    // check if something to update
    if (curCell === null) {
        return;
    }

    // add scores of ancestor cells
    for (var x = 0; x < curAncestor.parents.length; x++) {
        if(curAncestor.parents[x][0] > 0 && curAncestor.parents[x][1] > 0
            && curAncestor.parents[x][0] < curAncestor.parents[x][1]) {
            val *= this.getValue(curAncestor.parents[x][0], curAncestor.parents[x][1]) * Math.exp(1);
        }
    }
    // check for case (i,j-1)
    if(curAncestor.bps.length === 0){
        curCell.value = val;
    }
    else {
        curCell.value += val;
    }
    // curCell.traces.push(curAncestor);
}
;


/**
 * Fills the matrix according to the recursion.
 *
 * @param {string} sequence the RNA sequence to compute the matrix for
 * @param {int} minLoopLength the minimal loop length to be used for computation
 *
 * @returns {McKaskill_simple} this for call chaining
 */
McKaskill_simple.computeMatrix = function (sequence, minLoopLength) {

// resize and initialize matrix
    this.init(sequence, "structuresEnergy");
    //console.log(sequence, this.sequence.length);
// store minimal loop length
    this.minLoopLength = minLoopLength;

// fill matrix by diagonals
// iterate over all substructure spans that can have a base pair
    for (var span = minLoopLength; span < this.getDim(); span++) {
        // iterate over all rows
        for (var i = 0; i < this.getDim() - minLoopLength; i++) {
            // get column for current span
            var j = i + span;

            // j unpaired
            this.updateCell(i, j, Object.create(NussinovCellTrace).init([[i, j - 1]], []));

            // check base pair based decomposition : (k,j) base pair
            for (var k = i; k + minLoopLength < j; k++) {
                // check if sequence positions are compatible
                if (RnaUtil.areComplementary(this.sequence[k - 1], this.sequence[j - 1])) {
                    this.updateCell(i, j, Object.create(NussinovCellTrace).init([[i, k - 1], [k + 1, j - 1]], [[k, j]]));
                }
                ;
            }
            ;
        }
        ;
    }
    ;

    return this;
};



/****** McKaskill_base extending NussinovMatrix ************************/



/**
 * Implements a non-ambiguous recursion
 */
var McKaskill_base = Object.create(NussinovMatrix);

/**
 * Returns a description for the implemented recursion
 *
 * @returns {string} description of the recursion
 */
McKaskill_base.getRecursionDescription = function () {
    return "Recursion to count the energy of all the structures.";
};

/**
 * Access to the recursion in LaTeX encoding that is used in the computeMatrix implementation
 *
 * @returns {string} latex encoding of the recursion
 */
McKaskill_base.getRecursionInLatex = function () {
    return "$$Q_{i,j}^{b} = \\begin{cases} Q_{i + 1, j - 1} * \\exp(-E(bp)/RT) & \\text{ if }i,j \\text{ can form base pair} \\\\ 0 & \\text{ otherwise}\\end{cases}$$";
};

/**********  OVERWRITING + EXTENSION  ********************/

McKaskill_base.updateCell = function(i, j){
    // get cell to update
    var curCell = this.getCell(i, j);
    // check if something to update
    if (curCell === null) {
        return;
    }
    curCell.value = McKaskill_simple.getValue(i + 1, j - 1) * Math.exp(1);
}
;


/**
 * Fills the matrix according to the recursion.
 *
 * @param {string} sequence the RNA sequence to compute the matrix for
 * @param {int} minLoopLength the minimal loop length to be used for computation
 *
 * @returns {McKaskill_base} this for call chaining
 */
McKaskill_base.computeMatrix = function (sequence, minLoopLength) {

// resize and initialize matrix
    this.init(sequence, "structuresEnergy");
    //console.log(sequence, this.sequence.length);
// store minimal loop length
    this.minLoopLength = minLoopLength;

// fill matrix by diagonals
// iterate over all substructure spans that can have a base pair
    for (var span = minLoopLength; span < this.getDim(); span++) {
        // iterate over all rows
        for (var i = 0; i < this.getDim() - minLoopLength; i++) {
            // get column for current span
            var j = i + span;

            // j unpaired
            this.updateCell(i, j);

            // check base pair based decomposition : (k,j) base pair
            //for (var k = i; k + minLoopLength < j; k++) {
                // check if sequence positions are compatible
                //if (RnaUtil.areComplementary(this.sequence[k - 1], this.sequence[j - 1])) {
                    //this.updateCell(i, j, Object.create(NussinovCellTrace).init([[i, k - 1], [k + 1, j - 1]], [[k, j]]));
                //}
                ;
           // }
            ;

        }
        ;
    }
    ;

    return this;
};



/**
 * global list of available Nussinov algorithm implementations
 */
var availableAlgorithms = {

    /** ambiguous recursion */
    nussinovOriginal: NussinovMatrix_ambiguous,

    /** original unique recursion */
    nussinovUnique: NussinovMatrix_unique,

    /** Counting */
    counting: NussinovMatrix_structuresCount,

    /** McCaskill */
    mcKaskill: NussinovDPAlgorithm_McKaskill,

    /** counting3 */
    mycounting: NussinovDPAlgorithm_structuresCount,

};
//var availableAlgorithms = [NussinovMatrix_ambiguous, NussinovMatrix_unique];

/***************   TESTING  ********************/
//console.log(NussinovMatrix_structuresCount.getRecursionDescription());
//var m = NussinovMatrix_structuresCount.computeMatrix("CGGGC", 0);
//console.log(m.toString());
//
//{
//    var selectedAlgorithm = "nussinovOriginal";
////console.log("\n####  using", selectedAlgorithm, "  ####\n");
////console.log(availableAlgorithms[selectedAlgorithm].getRecursionDescription() + "\n");
////console.log(availableAlgorithms[selectedAlgorithm].getRecursionInLatex() + "\n");
////console.log(availableAlgorithms[selectedAlgorithm].computeMatrix("CGAGC", 0).toString());
////console.log(availableAlgorithms[selectedAlgorithm].computeMatrix("CGAGC", 1).toString());
////console.log("optimal structure :", availableAlgorithms[selectedAlgorithm].getOptimalTraceback().structure);
//
//    //selectedAlgorithm = "nussinovUnique";
//    console.log("\n####  using", selectedAlgorithm, "  ####\n");
//    //console.log(availableAlgorithms[selectedAlgorithm].getRecursionDescription() + "\n");
//    //console.log(availableAlgorithms[selectedAlgorithm].getRecursionInLatex() + "\n");
//    //console.log(availableAlgorithms[selectedAlgorithm].computeMatrix("GCGCG", 0).toString());
//    //console.log(availableAlgorithms[selectedAlgorithm].computeMatrix("GCGCG", 1).toString());
//    var xmat = availableAlgorithms[selectedAlgorithm].computeMatrix("CGAGCAGC", 0);
//    var opt = Object.create(Traceback);
//    opt.init(xmat.sequence.length);
//
//    var xy = wuchty_2nd(xmat, 0, selectedAlgorithm);
//    console.log("wuchty trace\n");
//
//    for (var i in xy) {
//        console.log(JSON.stringify(xy[i].traces), JSON.stringify(xy[i].structure));
//    }
//
//
//    //var x = xmat.traceallOpt(1, xmat.sequence.length);
//    //var struct = "";
//    //for (var i in x) {
//    //    struct += x[i].structure + ",";
//    //}
//    //console.log();
//    //for (var i in x) {
//        //console.log(JSON.stringify(x[i].traces), JSON.stringify(x[i].structure))
//    //}
//
//}
