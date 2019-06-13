##
# mex function builder.

class MexBuilder:
	
	# Build all mex functions
	@staticmethod
	def all():
		MexBuilder.solve_cohort()
		MexBuilder.generate_distribution()
		
	# Build solve_cohort
	@staticmethod
	def solve_cohort():
		MexBuilder.build('solve_cohort')
		
	# Build generate_distribution
	@staticmethod
	def generate_distribution():
		MexBuilder.build('generate_distribution')
		
	# Build mex function according to source function name
	@staticmethod
	def build(fname):
		
		print('\nBuilding mex function for %s.' % fname)
		
		# Specify code generation directory
		codegen_dir = '%s_codegen' % fname
		
		# Configure code generation
		mex_cfg = coder.config('mex')
		
		mex_cfg.ExtrinsicCalls            = false
		mex_cfg.IntegrityChecks           = false
		mex_cfg.SaturateOnIntegerOverflow = false
		
		# Generate mex function
		codegen('-d', codegen_dir, '-config', 'mex_cfg', '-o', fname, fname)
		
		# Clean up code generation directory
		rmdir(codegen_dir, 's')
		
		print('\nBuild complete.\n')