/*

    Copyright (C) 2017 Stanford HIVDB team

    Sierra is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Sierra is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

package edu.stanford.hivdb.graphql;

import graphql.schema.*;
import graphql.schema.GraphQLCodeRegistry.Builder;

import static graphql.Scalars.*;
import static graphql.schema.GraphQLObjectType.newObject;
import static graphql.schema.GraphQLCodeRegistry.newCodeRegistry;
import static graphql.schema.FieldCoordinates.coordinates;

import java.util.Collection;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.apache.commons.lang3.tuple.Triple;

import com.google.common.collect.Sets;

import edu.stanford.hivdb.drugresistance.GeneDR;
import edu.stanford.hivdb.drugresistance.algorithm.DrugResistanceAlgorithm;
import edu.stanford.hivdb.viruses.Gene;
import edu.stanford.hivdb.viruses.Virus;
import edu.stanford.hivdb.mutations.MutationSet;
import edu.stanford.hivdb.utilities.SimpleMemoizer;
import edu.stanford.hivdb.utilities.ValidationResult;

import static edu.stanford.hivdb.graphql.MutationSetDef.*;
import static edu.stanford.hivdb.graphql.DrugResistanceDef.*;
import static edu.stanford.hivdb.graphql.ValidationResultDef.*;
import static edu.stanford.hivdb.graphql.MutationPrevalenceDef.*;
import static edu.stanford.hivdb.graphql.AlgorithmComparisonDef.*;
import static edu.stanford.hivdb.graphql.DrugResistanceAlgorithmDef.*;

public class MutationsAnalysisDef {
	
	private static <VirusT extends Virus<VirusT>> Map<Gene<VirusT>, MutationSet<VirusT>> getMutationsByGeneFromSource(DataFetchingEnvironment env) {
		Triple<Set<Gene<VirusT>>, MutationSet<VirusT>, String> data = env.getSource();
		Set<Gene<VirusT>> knownGenes = data.getLeft();
		MutationSet<VirusT> mutations = data.getMiddle();
		
		Map<Gene<VirusT>, MutationSet<VirusT>> mutationsByGene = mutations.groupByGene();
		for (Gene<VirusT> gene : knownGenes) {
			if (!mutationsByGene.containsKey(gene)) {
				mutationsByGene.put(gene, new MutationSet<>());
			}
		}
		return mutationsByGene;
	}

	private static <VirusT extends Virus<VirusT>> MutationSet<VirusT> getMutationSetFromSource(DataFetchingEnvironment env) {
		Triple<Set<Gene<VirusT>>, MutationSet<VirusT>, String> data = env.getSource();
		return data.getMiddle();
	}
	
	private static <VirusT extends Virus<VirusT>> DataFetcher<List<ValidationResult>> makeMutValidationResultDataFetcher(VirusT virusIns) {
		return env -> {
			MutationSet<VirusT> mutations = getMutationSetFromSource(env);
			Collection<String> includeGenes = env.getArgument("includeGenes");
			return virusIns.validateMutations(mutations, Sets.newLinkedHashSet(includeGenes));
		};
	};

	private static DataFetcher<String> mutsNameDataFetcher = env -> {
		Triple<Set<Gene<?>>, MutationSet<?>, String> data = env.getSource();
		return data.getRight();
	};
	
	private static <VirusT extends Virus<VirusT>> DataFetcher<List<Map<String, Object>>> makeMutAllGeneMutSetDataFetcher(VirusT virusIns) {
		return env -> {
			Collection<String> includeGenes = env.getArgument("includeGenes");
			Set<String> includeGeneSet = Sets.newHashSet(includeGenes);
			return getMutationsByGeneFromSource(env)
				.entrySet()
				.stream()
				.map(entry -> Map.of(
					"gene", (Object) entry.getKey(),
					"mutations", (Object) entry.getValue()
				))
				.filter(map -> includeGeneSet.contains(((Gene<?>) map.get("gene")).getAbstractGene()))
				.collect(Collectors.toList());
		};
	}

	private static <VirusT extends Virus<VirusT>> DataFetcher<List<GeneDR<VirusT>>> makeMutDRDataFetcher(VirusT virusIns) {
		return env -> {
			Map<Gene<VirusT>, MutationSet<VirusT>> mutationsByGene = getMutationsByGeneFromSource(env);
			String algName = env.getArgument("algorithm");
			Collection<String> includeGenes = env.getArgument("includeGenes");
			Set<String> includeGeneSet = Sets.newHashSet(includeGenes);
			return mutationsByGene
				.entrySet()
				.stream()
				.map(e -> new GeneDR<>(
					e.getKey(), e.getValue(), virusIns.getDrugResistAlgorithm(algName)
				))
				.filter(geneDR -> includeGeneSet.contains(geneDR.getAbstractGene()))
				.collect(Collectors.toList());
		};
	};
	
	private static DataFetcher<List<Map<String, Object>>> boundMutPrevListDataFetcher = env -> {
		MutationSet<?> mutations = getMutationSetFromSource(env);
		Collection<String> includeGenes = env.getArgument("includeGenes");
		return getBoundMutationPrevalenceList(mutations, Sets.newHashSet(includeGenes));
	};
	
	private static <VirusT extends Virus<VirusT>> DataFetcher<List<Map<String, Object>>> makeMutAlgCmpDataFetcher(VirusT virusIns) {
		return env -> {
			List<String> asiAlgs = env.getArgument("algorithms");
			List<Map<String, String>> customAlgs = env.getArgument("customAlgorithms");
			if (asiAlgs == null) { asiAlgs = Collections.emptyList(); }
			if (customAlgs == null) { customAlgs = Collections.emptyList(); }
			if (asiAlgs.isEmpty() && customAlgs.isEmpty()) {
				return Collections.emptyList();
			}
			asiAlgs = asiAlgs
				.stream().filter(alg -> alg != null)
				.collect(Collectors.toList());
			Map<String, String> customAlgs2 = customAlgs
				.stream()
				.filter(map -> map != null)
				.collect(Collectors.toMap(
					map -> map.get("name"),
					map -> map.get("xml"),
					(x1, x2) -> x2,
					LinkedHashMap::new
				));
			MutationSet<VirusT> mutations = getMutationSetFromSource(env);
			return fetchAlgorithmComparisonData(virusIns, mutations, asiAlgs, customAlgs2);
		};
	};

	public static <VirusT extends Virus<VirusT>> GraphQLCodeRegistry makeMutationsAnalysisCodeRegistry(VirusT virusIns) {
		Builder builder = newCodeRegistry()
			.dataFetcher(
				coordinates("MutationsAnalysis", "name"),
				mutsNameDataFetcher
			)
			.dataFetcher(
				coordinates("MutationsAnalysis", "validationResults"),
				makeMutValidationResultDataFetcher(virusIns)
			)
			.dataFetcher(
				coordinates("MutationsAnalysis", "allGeneMutations"),
				makeMutAllGeneMutSetDataFetcher(virusIns)
			)
			.dataFetcher(
				coordinates("MutationsAnalysis", "mutationPrevalences"),
				boundMutPrevListDataFetcher
			);
		
		DrugResistanceAlgorithm<?> defaultDRAlgo = virusIns.getDefaultDrugResistAlgorithm();
		if (defaultDRAlgo != null) {
			builder = builder
				.dataFetcher(
					coordinates("MutationsAnalysis", "drugResistance"),
					makeMutDRDataFetcher(virusIns)
				)
				.dataFetcher(
					coordinates("MutationsAnalysis", "algorithmComparison"),
					makeMutAlgCmpDataFetcher(virusIns)
				);
		
		}
		return builder.build();
	};

	public static SimpleMemoizer<GraphQLObjectType> oGeneMutations = new SimpleMemoizer<>(
		virusName -> newObject()
			.name("GeneMutations")
			.field(field -> field
				.type(new GraphQLTypeReference("Gene"))
				.name("gene")
				.description("Gene of the mutation set.")
			)
			.field(field -> newMutationSet(virusName, field, "mutations", /* enableIncludedGenes= */true)
				.description("All mutations of this gene.")
			)
			.build()
	);
	
	public static SimpleMemoizer<GraphQLObjectType> oMutationsAnalysis = new SimpleMemoizer<>(
		virusName -> {
			GraphQLObjectType.Builder builder = newObject()
			.name("MutationsAnalysis")
			.field(field -> field
				.type(GraphQLString)
				.name("name")
				.description("Optional name provided by client to identify this mutation list."))
			.field(field -> field
				.type(new GraphQLList(oValidationResult))
				.name("validationResults")
				.argument(arg -> arg
					.type(new GraphQLList(GeneDef.enumGene.get(virusName)))
					.name("includeGenes")
					.defaultValue(Virus.getInstance(virusName).getDefaultIncludedGenes())
					.description("Genes to be included in the results")
				)
				.description("Validation results for the mutation list."))
			.field(field -> field
				.type(new GraphQLList(oGeneMutations.get(virusName)))
				.name("allGeneMutations")
				.argument(arg -> arg
					.type(new GraphQLList(GeneDef.enumGene.get(virusName)))
					.name("includeGenes")
					.defaultValue(Virus.getInstance(virusName).getDefaultIncludedGenes())
					.description("Genes to be included in the results")
				)
				.description("Mutations groupped by gene."))
			.field(field -> field
				.type(new GraphQLList(oBoundMutationPrevalence.get(virusName)))
				.name("mutationPrevalences")
				.argument(arg -> arg
					.type(new GraphQLList(GeneDef.enumGene.get(virusName)))
					.name("includeGenes")
					.defaultValue(Virus.getInstance(virusName).getDefaultIncludedGenes())
					.description("Genes to be included in the results")
				)
				.description("List of mutation prevalence results."));
			Virus<?> virusIns = Virus.getInstance(virusName);
			DrugResistanceAlgorithm<?> defaultDRAlgo = virusIns.getDefaultDrugResistAlgorithm();
			if (defaultDRAlgo != null) {
				builder = builder
					.field(field -> field
						.type(new GraphQLList(oDrugResistance.get(virusName)))
						.name("drugResistance")
						.argument(arg -> arg
							.name("algorithm")
							.type(oASIAlgorithm.get(virusName))
							.defaultValue(Virus.getInstance(virusName).getDefaultDrugResistAlgorithm().getName())
							.description("One of the built-in ASI algorithms."))
						.argument(arg -> arg
							.type(new GraphQLList(GeneDef.enumGene.get(virusName)))
							.name("includeGenes")
							.defaultValue(Virus.getInstance(virusName).getDefaultIncludedGenes())
							.description("Genes to be included in the results")
						)
						.description("List of drug resistance results by genes."))
					.field(field -> field
						.type(new GraphQLList(oAlgorithmComparison.get(virusName)))
						.name("algorithmComparison")
						.description("List of ASI comparison results.")
						.argument(aASIAlgorithmArgument.get(virusName))
						.argument(aASICustomAlgorithmArgument));
			}
			builder = virusIns.getVirusGraphQLExtension().extendObjectBuilder("MutationsAnalysis", builder);
			return builder.build();
		}
	);

}
