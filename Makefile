all:
	@echo "> Compiling Simulation Results.\n"
	+$(MAKE) -C Coding/"Python Simulations"
	@echo "> Moving results to dissertation project.\n"
	+$(MAKE) -C Results

qiskit:
	@echo "> Compiling Qiskit Results.\n"
	+$(MAKE) -C Coding/"Qiskit"
	@echo "> Moving results to dissertation project.\n"
	+$(MAKE) -C Results
git:
	git branch
