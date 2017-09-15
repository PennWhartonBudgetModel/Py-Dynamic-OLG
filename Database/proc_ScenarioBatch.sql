/**
	p_ScenarioBatch

	Read the list of DynamicModel Scenarios with XML fields
	pulled out into columns.

**/
USE second_chart;

SET NOCOUNT ON;


IF OBJECT_ID(N'p_ScenarioBatch', N'P' ) IS NOT NULL
	DROP PROCEDURE p_ScenarioBatch;

GO
;
CREATE PROCEDURE p_ScenarioBatch
	@BatchID int
AS
BEGIN
	SELECT	BatchID
		,	ID
		,	OpenEconomy
		,	DynamicSIMParameters.value('(/Dynamics/SavingsElasticity)[1]', 'real') AS SavingsElasticity
		,	DynamicSIMParameters.value('(/Dynamics/LaborElasticity)[1]', 'real') AS LaborElasticity
		,	DynamicSIMParameters.value('(/Dynamics/ExpenditureShift)[1]', 'real') AS ExpenditureShift
		,	DynamicSIMParameters.value('(/Dynamics/UseDynamicBaseline)[1]', 'bit') AS UseDynamicBaseline
		,	DynamicSIMParameters.value('(/Dynamics/TaxPlan)[1]', 'varchar(20)') AS TaxPlan
	FROM Scenario
	WHERE BatchID = @BatchID
		AND ModelType = 'D'
	;
END

GO
;

EXEC u_SetProcPermissions

PRINT 'Built p_ScenarioBatch';

/* END -- script */