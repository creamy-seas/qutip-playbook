OSFLAG :=
ifeq ($(OS),Windows_NT)
	OSFLAG += WIN32
else
	OSFLAG += UNIX
endif

jupyter:
	@docker-compose up --build

ilya-testing: $(OSFLAG)
	@echo "Do not use - this is WIP for easier deployment"

UNIX:
	@echo "⚛  Assuming UNIX-like environment"
	@echo "⚙ Building system - this may take 5 minutes"
	@docker build -t rhul-qutip-notebook .
	@echo docker run --publish 8888:8888\
						--name rhul-qutip-notebook-container\
						--rm\
						-i\
						-t\
						--volume "$pwd":/qutip-simulator rhul-qutip-notebook

WIN:
	@echo "⚛  Assuming WINDOWS-like environment"
	@echo "⚙ Building system - this may take 5 minutes"
	@docker build -t rhul-qutip-notebook .
	@docker run --publish 8888:8888\
						--name rhul-qutip-notebook-container\
						--rm\
						-i\
						-t\
						--volume ${PWD}:/qutip-simulator\
						rhul-qutip-notebook
