project(alpha_beta)

file(GLOB ALPHA_BETA_SOURCE
	"*.hpp"
	"*.cpp"
	)

add_executable(alpha_beta ${ALPHA_BETA_SOURCE})
target_link_libraries(alpha_beta saedb)
