#include <vector>
#include <iosfwd>

namespace gr {
	namespace trellis {

	// Finite State Machine Specification class (trellis_coding_blk): An instance of this class represents a finite state machine specification (FSMS)
	// rather than the FSM itself. In particular the state of the FSM is not stored within an instance of this class.

	class fsm
	{
	private:
		int d_I;				// input alphabet cardinality.
		int d_S;				// number of states.
		int d_O;				// output alphabet cardinality.
		std::vector<int> d_NS;			// next_state (NS) = d_NS[current_state * d_I + input_symbol]
		std::vector<int> d_OS;			// output_symbol (OS) = d_OS[current_state * d_I + input_symbol]
		std::vector< std::vector<int> > d_PS;	// previous state (PS)

		// d_PS[current_state][k] & d_PI[current_state][k], are pairs of form (prev_state, prev_input_symbol) that could produce the current state.
		std::vector< std::vector<int> > d_PI; // previous input symbol (PI)

		// d_TMl[s*d_S+es] is the shortest number of steps to get from state s to state es.
		std::vector<int> d_TMl; // termination matrix (TM)

		// d_TMi[s*d_S+es] is the input symbol required to set off on the shortest path from state s to es.
		std::vector<int> d_TMi;

		void generate_PS_PI ();
		void generate_TM ();
		bool find_es(int es);

	public:
		fsm(); //brief Constructor to create an uninitialized FSMS.
		fsm(const fsm &FSM); //brief Constructor to copy an FSMS.

		/*!
		* \brief Constructor to to create an FSMS.
		*
		* \param I	     The number of possible input symbols.
		* \param S           The number of possible FSM states.
		* \param O           The number of possible output symbols.
		* \param NS          A mapping from (current state, input symbol) to next state.
		*                    next_state = NS[current_state * I + input_symbol]
		* \param OS          A mapping from (current state, input symbol) to output symbol.
		*                    output_symbol = OS[current_state * I + input_symbol]
		*
		*/
		fsm(int I, int S, int O, const std::vector<int> &NS, const std::vector<int> &OS);

		// brief Constructor to create an FSMS from file contents.
		fsm(const char *name); // uses a file as input, enter filename



		// brief Creates an FSMS from the generator matrix of a (n, k) binary convolutional code.
		fsm(int k, int n, const std::vector<int> &G); // <-- THIS ONE



		/*!
		* \brief Creates an FSMS describing ISI.
		*
		* \param mod_size    modulation size
		* \param ch_length   channel length
		*
		*/
		fsm(int mod_size, int ch_length);

		/*!
		* \brief Creates an FSMS describing the trellis for a CPM.
		*
		* \param P    ???? h=K/P (relatively prime)
		* \param M    alphabet size
		* \param L    pulse duration
		*
		* This FSM is based on the paper by B. Rimoldi
		* "A decomposition approach to CPM", IEEE Trans. Info Theory, March 1988
		* See also my own notes at http://www.eecs.umich.edu/~anastas/docs/cpm.pdf
		*/
		fsm(int P, int M, int L);

		/*!
		* \brief Creates an FSMS describing the joint trellis of two FSMs.
		*
		* \param FSM1  first FSMS
		* \param FSM2  second FSMS
		*/
		fsm(const fsm &FSM1, const fsm &FSM2);

		/*!
		* \brief Creates an FSMS representing n stages through the originial FSM (AKA radix-n FSM).
		*
		* \param FSM      Original FSMs
		* \param n        Number of stages.
		*/
		fsm(const fsm &FSM, int n);
		int I() const { return d_I; }
		int S() const { return d_S; }
		int O() const { return d_O; }
		const std::vector<int> & NS() const { return d_NS; }
		const std::vector<int> & OS() const { return d_OS; }
		const std::vector< std::vector<int> > & PS() const { return d_PS; }
		const std::vector< std::vector<int> > & PI() const { return d_PI; }
		const std::vector<int> & TMi() const { return d_TMi; }
		const std::vector<int> & TMl() const { return d_TMl; }

		/*!
		* \brief Creates an svg image of the trellis representation.
		*
		* \param filename         filename
		* \param number_stages    ????
		*/
		void write_trellis_svg(std::string filename ,int number_stages);

		/*!
		* \brief Write the FSMS to a file.
		*
		* \param filename   filename
		*/
		void write_fsm_txt(std::string filename);
	};

	} /* namespace trellis */
} /* namespace gr */
