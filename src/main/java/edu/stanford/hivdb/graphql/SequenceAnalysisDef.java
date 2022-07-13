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
import static graphql.Scalars.*;
import static graphql.schema.GraphQLArgument.newArgument;
import static graphql.schema.GraphQLObjectType.newObject;
import static graphql.schema.GraphQLCodeRegistry.newCodeRegistry;
import static graphql.schema.FieldCoordinates.coordinates;

import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;

import edu.stanford.hivdb.drugresistance.GeneDR;
import edu.stanford.hivdb.genotypes.BoundGenotype;
import edu.stanford.hivdb.genotypes.GenotypeResult;
import edu.stanford.hivdb.mutations.FrameShift;
import edu.stanford.hivdb.mutations.MutationSet;
import edu.stanford.hivdb.sequences.AlignedGeneSeq;
import edu.stanford.hivdb.sequences.AlignedSequence;
import edu.stanford.hivdb.utilities.SimpleMemoizer;
import edu.stanford.hivdb.utilities.ValidationResult;
import edu.stanford.hivdb.viruses.Virus;

import static edu.stanford.hivdb.graphql.UnalignedSequenceDef.*;
import static edu.stanford.hivdb.graphql.MutationSetDef.*;
import static edu.stanford.hivdb.graphql.GeneDef.*;
import static edu.stanford.hivdb.graphql.StrainDef.*;
import static edu.stanford.hivdb.graphql.FrameShiftDef.*;
import static edu.stanford.hivdb.graphql.SubtypeDef.*;
import static edu.stanford.hivdb.graphql.AlignedGeneSequenceDef.*;
import static edu.stanford.hivdb.graphql.DrugResistanceDef.*;
import static edu.stanford.hivdb.graphql.ValidationResultDef.*;
import static edu.stanford.hivdb.graphql.SubtypeV2Def.*;
import static edu.stanford.hivdb.graphql.MutationPrevalenceDef.*;
import static edu.stanford.hivdb.graphql.AlgorithmComparisonDef.*;
import static edu.stanford.hivdb.graphql.DrugResistanceAlgorithmDef.*;

public class SequenceAnalysisDef {
	
	private static <VirusT extends Virus<VirusT>> DataFetcher<List<Map<String, Object>>> makeSubtypesDataFetcher(VirusT virusIns) {
		return env -> {
			int first = env.getArgument("first");
			AlignedSequence<VirusT> alignedSeq = env.getSource();
			GenotypeResult<VirusT> subtypeResult = alignedSeq.getGenotypeResult();
			if (subtypeResult == null) {
				return Collections.emptyList();
			}
			List<BoundGenotype<VirusT>> matches = subtypeResult.getAllMatches();
			first = Math.min(first, matches.size());
			List<BoundGenotype<VirusT>> subtypes = matches.subList(0, first);
			return subtypes
			.stream()
			.map(g -> {
				Map<String, Object> r = new HashMap<>();
				String distancePcnt = g.getDistancePcnt();
				distancePcnt = distancePcnt.substring(0, distancePcnt.length() - 1);
				r.put("name", g.getGenotype().getIndexName());
				r.put("distancePcnt", Double.parseDouble(distancePcnt));
				r.put("display", g.getDisplay());
				return r;
			})
			.collect(Collectors.toList());
		};
	};

	private static <VirusT extends Virus<VirusT>> DataFetcher<List<BoundGenotype<VirusT>>> makeSubtypesDataFetcherV2(VirusT virusIns) {
		return env -> {
			int first = env.getArgument("first");
			AlignedSequence<VirusT> alignedSeq = env.getSource();
			GenotypeResult<VirusT> subtypeResult = alignedSeq.getGenotypeResult();
			if (subtypeResult == null) {
				return Collections.emptyList();
			}
			List<BoundGenotype<VirusT>> matches = subtypeResult.getAllMatches();
			first = Math.min(first, matches.size());
			return matches.subList(0, first);
		};
	};

	private static <VirusT extends Virus<VirusT>> DataFetcher<List<ValidationResult>> makeValidationResultsDataFetcher(VirusT virusIns) {
		return env -> {
			AlignedSequence<VirusT> alignedSeq = env.getSource();
			Collection<String> includeGenes = env.getArgument("includeGenes");
			return alignedSeq.getValidationResults(Sets.newLinkedHashSet(includeGenes));
		};
	}
	
	private static <VirusT extends Virus<VirusT>> DataFetcher<List<AlignedGeneSeq<VirusT>>> makeAlignedGeneSequencesDataFetcher(VirusT virusIns) {
		return env -> {
			AlignedSequence<VirusT> alignedSeq = env.getSource();
			Collection<String> includeGenes = env.getArgument("includeGenes");
			return alignedSeq.getAlignedGeneSequences(Sets.newLinkedHashSet(includeGenes));
		};
	}

	private static <VirusT extends Virus<VirusT>> DataFetcher<List<GeneDR<VirusT>>> makeDrugResistanceDataFetcher(VirusT virusIns) {
		return env -> {
			AlignedSequence<VirusT> alignedSeq = env.getSource();
			String algName = env.getArgument("algorithm");
			Collection<String> includeGenes = env.getArgument("includeGenes");
			List<AlignedGeneSeq<VirusT>> geneSeqs = alignedSeq.getAlignedGeneSequences(Sets.newLinkedHashSet(includeGenes));
			return Lists.newArrayList(
				GeneDR.newFromAlignedGeneSeqs(
					geneSeqs, virusIns.getDrugResistAlgorithm(algName)
				).values()
			);
		};
	};
	
	private static DataFetcher<List<Map<String, Object>>> boundMutPrevListDataFetcher = env -> {
		AlignedSequence<?> alignedSeq = env.getSource();
		MutationSet<?> mutations = alignedSeq.getMutations();
		Collection<String> includeGenes = env.getArgument("includeGenes");
		return getBoundMutationPrevalenceList(mutations, Sets.newHashSet(includeGenes));
	};

	private static <VirusT extends Virus<VirusT>> DataFetcher<List<FrameShift<VirusT>>> makeFrameShiftsDataFetcher(VirusT virusIns) {
		return env -> {
			AlignedSequence<VirusT> alignedSeq = env.getSource();
			List<FrameShift<VirusT>> fss = alignedSeq.getFrameShifts();
			Collection<String> includeGenes = env.getArgument("includeGenes");
			Set<String> includeGeneSet = Sets.newHashSet(includeGenes);
			return fss.stream()
				.filter(fs -> includeGeneSet.contains(fs.getAbstractGene()))
				.collect(Collectors.toList());
		};
	};
	
	private static <VirusT extends Virus<VirusT>> DataFetcher<List<Map<String, Object>>> makeAlgComparisonDataFetcher(VirusT virusIns) {
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
			AlignedSequence<VirusT> alignedSeq = env.getSource();
			return fetchAlgorithmComparisonData(virusIns, alignedSeq.getMutations(), asiAlgs, customAlgs2);
		};
	};
			
	public static <VirusT extends Virus<VirusT>> GraphQLCodeRegistry makeSequenceAnalysisCodeRegistry(VirusT virusIns) {
		return newCodeRegistry()
			.dataFetcher(
				coordinates("SequenceAnalysis", "subtypes"),
				makeSubtypesDataFetcher(virusIns)
			)
			.dataFetcher(
				coordinates("SequenceAnalysis", "subtypesV2"),
				makeSubtypesDataFetcherV2(virusIns)
			)
			.dataFetcher(
				coordinates("SequenceAnalysis", "genotypes"),
				makeSubtypesDataFetcherV2(virusIns)
			)
			.dataFetcher(
				coordinates("SequenceAnalysis", "validationResults"),
				makeValidationResultsDataFetcher(virusIns)
			)
			.dataFetcher(
				coordinates("SequenceAnalysis", "alignedGeneSequences"),
				makeAlignedGeneSequencesDataFetcher(virusIns)
			)
			.dataFetcher(
				coordinates("SequenceAnalysis", "drugResistance"),
				makeDrugResistanceDataFetcher(virusIns)
			)
			.dataFetcher(
				coordinates("SequenceAnalysis", "mutationPrevalences"),
				boundMutPrevListDataFetcher
			)
			.dataFetcher(
				coordinates("SequenceAnalysis", "algorithmComparison"),
				makeAlgComparisonDataFetcher(virusIns)
			)
			.dataFetcher(
				coordinates("SequenceAnalysis", "mutations"),
				newMutationSetDataFetcher(virusIns, "mutations")
			)
			.dataFetcher(
				coordinates("SequenceAnalysis", "frameShifts"),
				makeFrameShiftsDataFetcher(virusIns)
			)
			.dataFetchers(makeAlignedGeneSequenceCodeRegistry(virusIns))
			.build();
	}

	public static SimpleMemoizer<GraphQLObjectType> oSequenceAnalysis = new SimpleMemoizer<>(
		virusName -> {
			GraphQLObjectType.Builder builder = newObject()
				.name("SequenceAnalysis")
					.field(field -> field
						.type(oUnalignedSequence)
						.name("inputSequence")
						.description("The original unaligned sequence."))
					.field(field -> field
						.type(oStrain)
						.name("strain")
						.description("Virus strain of this sequence."))
					.field(field -> field
						.type(GraphQLBoolean)
						.name("isReverseComplement")
						.description("True if the alignment result was based on the reverse complement of input sequence."))
					.field(field -> field
						.type(new GraphQLList(oGene.get(virusName)))
						.name("availableGenes")
						.description("Available genes found in the sequence."))
					.field(field -> field
						.type(new GraphQLList(oValidationResult))
						.name("validationResults")
						.argument(arg -> arg
							.type(new GraphQLList(GeneDef.enumGene.get(virusName)))
							.name("includeGenes")
							.defaultValue(Virus.getInstance(virusName).getDefaultIncludedGenes())
							.description("Genes to be included in the results")
						)
						.description("Validation results for this sequence."))
					.field(field -> field
						.type(new GraphQLList(oAlignedGeneSequence.get(virusName)))
						.name("alignedGeneSequences")
						.argument(arg -> arg
							.type(new GraphQLList(GeneDef.enumGene.get(virusName)))
							.name("includeGenes")
							.defaultValue(Virus.getInstance(virusName).getDefaultIncludedGenes())
							.description("Genes to be included in the results")
						)
						.description("List of aligned sequence distinguished by genes."))
					.field(field -> field
						.type(new GraphQLList(oBoundSubtypeV2))
						.name("subtypesV2")
						.argument(newArgument()
							.type(GraphQLInt)
							.name("first")
							.defaultValue(2)
							.description(
								"Fetch only the first nth closest subtypes. Default to 2.")
							.build())
						.description(
							"List of virus groups or subtypes, or species. " +
							"Sorted by the similarity from most to least."))
					.field(field -> field
						.type(oBoundSubtypeV2)
						.name("bestMatchingSubtype")
						.description(
							"The best matching subtype."))
					.field(field -> field
						.type(new GraphQLList(oBoundSubtypeV2))
						.name("genotypes")
						.argument(arg -> arg
							.type(GraphQLInt)
							.name("first")
							.defaultValue(2)
							.description(
								"Fetch only the first nth closest genotypes. Default to 2.")
						)
						.deprecate("Use field `subtypesV2` instead.")
						.description(
							"List of virus groups or subtypes, or species. " +
							"Sorted by the similarity from most to least."))
					.field(field -> field
						.type(oBoundSubtypeV2)
						.name("bestMatchingGenotype")
						.deprecate("Use field `bestMatchingSubtype` instead.")
						.description(
							"The best matching genotype."))
					.field(field -> field
						.type(GraphQLFloat)
						.name("mixturePcnt")
						.deprecate("Use field `mixtureRate` * 100 instead."))
					.field(field -> field
						.type(GraphQLFloat)
						.name("mixtureRate")
						.description(
							"Mixture rate of the sequence. Notes only RYMWKS " +
							"are counted."))
					.field(field -> newMutationSet(virusName, field, "mutations", /* enableIncludedGenes= */true)
						.description("All mutations found in the aligned sequence."))

					.field(field -> field
						.type(GraphQLInt)
						.name("mutationCount")
						.description("Number of mutations without counting unsequenced regions and multiple continuous deletions"))
					.field(field -> field
						.type(GraphQLInt)
						.name("unusualMutationCount")
						.description("Number of unusual mutations without counting unsequenced regions and multiple continuous deletions"))
					.field(field -> field
						.type(GraphQLInt)
						.name("insertionCount")
						.description("Number of insertions"))
					.field(field -> field
						.type(GraphQLInt)
						.name("deletionCount")
						.description("Number of deletions"))
					.field(field -> field
						.type(GraphQLInt)
						.name("stopCodonCount")
						.description("Number of positions with stop codons without counting unsequenced regions and multiple continuous deletions"))
					.field(field -> field
						.type(GraphQLInt)
						.name("ambiguousMutationCount")
						.description("Number of ambiguous positions without counting unsequenced regions and multiple continuous deletions"))
					.field(field -> field
						.type(GraphQLInt)
						.name("apobecMutationCount")
						.description("Number of APOBEC mutations without counting unsequenced regions and multiple continuous deletions"))
					.field(field -> field
						.type(GraphQLInt)
						.name("apobecDRMCount")
						.description("Number of APOBEC DRMs without counting unsequenced regions and multiple continuous deletions"))
					
					.field(field -> field
						.type(GraphQLInt)
						.name("frameShiftCount")
						.description("Number of frame shifts"))
					.field(field -> field
						.type(new GraphQLList(oFrameShift.get(virusName)))
						.name("frameShifts")
						.argument(arg -> arg
							.type(new GraphQLList(GeneDef.enumGene.get(virusName)))
							.name("includeGenes")
							.defaultValue(Virus.getInstance(virusName).getDefaultIncludedGenes())
							.description("Genes to be included in the results")
						)
						.description("List of all frame shifts"))
					.field(field -> field
						.type(new GraphQLList(oDrugResistance.get(virusName)))
						.name("drugResistance")
						.argument(arg -> arg
							.name("algorithm")
							.type(oASIAlgorithm.get(virusName))
							.defaultValue(Virus.getInstance(virusName).getDefaultDrugResistAlgorithm().getName())
							.description("One of the built-in ASI algorithms.")
						)
						.argument(arg -> arg
							.type(new GraphQLList(GeneDef.enumGene.get(virusName)))
							.name("includeGenes")
							.defaultValue(Virus.getInstance(virusName).getDefaultIncludedGenes())
							.description("Genes to be included in the results")
						)
						.description("List of drug resistance results by genes."))
					.field(field -> field
						.type(new GraphQLList(oBoundMutationPrevalence.get(virusName)))
						.name("mutationPrevalences")
						.argument(arg -> arg
							.type(new GraphQLList(GeneDef.enumGene.get(virusName)))
							.name("includeGenes")
							.defaultValue(Virus.getInstance(virusName).getDefaultIncludedGenes())
							.description("Genes to be included in the results")
						)
						.description("List of mutation prevalence results."))
					.field(field -> field
						.type(new GraphQLList(oBoundSubtype.get(virusName)))
						.name("subtypes")
						.deprecate("Use field `subtypesV2` instead.")
						.argument(newArgument()
							.type(GraphQLInt)
							.name("first")
							.defaultValue(2)
							.description(
							"Fetch only the first nth closest subtypes. Default to 2.")
						.build())
					.description(
						"List of virus groups or subtypes, or species. " +
						"Sorted by the similarity from most to least."))
					.field(field -> field
						.type(GraphQLString)
						.name("subtypeText")
						.deprecate("Use field `bestMatchingSubtype { display }` instead.")
						.description(
							"Formatted text for best matching subtype."))
					.field(field -> field
						.type(new GraphQLList(oAlgorithmComparison.get(virusName)))
						.name("algorithmComparison")
						.description("List of ASI comparison results.")
						.argument(aASIAlgorithmArgument.get(virusName))
						.argument(aASICustomAlgorithmArgument));
			Virus<?> virusIns = Virus.getInstance(virusName);
			builder = virusIns.getVirusGraphQLExtension().extendObjectBuilder("SequenceAnalysis", builder);
			return builder.build();
		}
	);

}
