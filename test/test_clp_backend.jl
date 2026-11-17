using Test
using ReactiveModel
import ReactiveModel as RM

@testset "ClpBackend : rechargement et limite de temps" begin
    @model ClpSensProbe(data, sel) begin
        G = 1:2
        @variable p[g in G] in [0.0, 10.0]
        @constraint Demande begin
            sum(p[g] for g in G) == sel.demande
        end
        @objective Min begin
            sum(data.cout[g] * p[g] for g in G)
        end
    end

    donnees = (cout = [1.0, 3.0],)
    definition = ClpSensProbe(donnees)
    lp = RM.LinearProblem(RM.ConstData(definition, donnees), (demande = 12.0,); backend=ClpBackend())

    # Min : 10 unités au coût 1, puis 2 au coût 3.
    @test RM.objective_value(RM.solve(lp, (demande = 12.0,))) ≈ 16.0

    # Le rechargement ne doit pas transformer une minimisation en maximisation (ici 34.0).
    RM._invalidate_solution!(lp)
    lp.parameter_snapshot = RM.Parameters((demande = -1.0,))
    lp.workspace.last_status = :not_solved
    RM.backend_reset_basis!(lp.backend)
    @test RM.objective_value(RM.solve(lp, (demande = 12.0,))) ≈ 16.0

    # `Clp_setMaximumSeconds` est une échéance absolue en temps CPU du processus : ne jamais l'armer.
    @test ccall((:Clp_maximumSeconds, RM._CLP_LIB), Cdouble, (Ptr{Cvoid},), lp.backend.inner) < 0
end

@testset "ClpBackend : le rechargement garde les coefficients mis à jour" begin
    # Un coefficient de membre gauche qui dépend des paramètres : `a * v[1] + v[2] >= 4`.
    @model ClpCoefProbe(data, sel) begin
        I = 1:2
        @variable v[i in I] in [0.0, 10.0]
        @constraint Lien begin
            sel.a * v[1] + v[2] >= 4.0
        end
        @objective Min begin
            v[1] + 2.0 * v[2]
        end
    end

    donnees = (k = 1.0,)
    lp = RM.LinearProblem(RM.ConstData(ClpCoefProbe(donnees), donnees), (a = 1.0,); backend=ClpBackend())
    @test RM.objective_value(RM.solve(lp, (a = 1.0,))) ≈ 4.0     # v = (4, 0)
    @test RM.objective_value(RM.solve(lp, (a = 0.0,))) ≈ 8.0     # v[2] >= 4

    # Remise à zéro de la base entre deux résolutions. Le
    # coefficient de v[1] vaut déjà 0 dans la copie de ReactiveModel : il n'est pas renvoyé. Si le rechargement
    # reprenait la matrice d'`initialize!` (a = 1), Clp résoudrait v[1] + v[2] >= 4 et trouverait 4.
    RM._invalidate_solution!(lp)
    lp.parameter_snapshot = RM.Parameters((a = -1.0,))
    lp.workspace.last_status = :not_solved
    RM.backend_reset_basis!(lp.backend)
    @test RM.objective_value(RM.solve(lp, (a = 0.0,))) ≈ 8.0
end
