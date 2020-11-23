all:
	@echo "Executing random walk models and then copying them to the Latex project"
	+$(MAKE) -C Coding/"Python Simulations"
	+$(MAKE) -C Results

git:
	git branch
