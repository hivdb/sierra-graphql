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
import static graphql.schema.GraphQLObjectType.newObject;
import static graphql.schema.GraphQLCodeRegistry.newCodeRegistry;
import static graphql.schema.FieldCoordinates.coordinates;

import edu.stanford.hivdb.drugs.DrugClass;
import edu.stanford.hivdb.utilities.SimpleMemoizer;
import edu.stanford.hivdb.viruses.Virus;

import static edu.stanford.hivdb.graphql.MutationSetDef.*;
import static edu.stanford.hivdb.graphql.ExtGraphQL.*;

public class DrugClassDef {

	public static GraphQLCodeRegistry drugClassCodeRegistry = newCodeRegistry()
		.dataFetcher(
			coordinates("DrugClass", "hasDrugResistMutations"),
			new ExtPropertyDataFetcher<String>("hasDrugResistMutations")
		)
		.dataFetcher(
			coordinates("DrugClass", "hasSurveilDrugResistMutations"),
			new ExtPropertyDataFetcher<String>("hasSurveilDrugResistMutations")
		)
		.dataFetcher(
			coordinates("DrugClass", "hasRxSelectedMutations"),
			new ExtPropertyDataFetcher<String>("hasRxSelectedMutations")
		)
		.build();

	public static SimpleMemoizer<GraphQLEnumType> oDrugClassEnum = new SimpleMemoizer<>(
		name -> {
			Virus<?> virusIns = Virus.getInstance(name);
			GraphQLEnumType.Builder
				newDrugClassEnum = GraphQLEnumType.newEnum()
				.name("DrugClassEnum");
			for (DrugClass<?> drugClass : virusIns.getDrugClasses()) {
				String dcText = drugClass.toString();
				newDrugClassEnum.value(dcText, dcText);
			}
			return newDrugClassEnum.build();
		}
	);

	public static SimpleMemoizer<GraphQLObjectType> oDrugClass = new SimpleMemoizer<>(
		name -> (
			newObject()
			.name("DrugClass")
			.description("Antiviral drug class.")
			.field(field -> field
				.type(oDrugClassEnum.get(name))
				.name("name")
				.description("Name of the drug class."))
			.field(field -> field
				.type(GraphQLString)
				.name("fullName")
				.description("Full name of the drug class."))
			.field(field -> field
				.type(new GraphQLList(new GraphQLTypeReference("Drug")))
				.name("drugs")
				.description("Drugs of this drug class."))
			.field(field -> field
				.type(new GraphQLTypeReference("Gene"))
				.name("gene")
				.description("Gene the drug class belongs to."))
			.field(field -> newMutationSet(name, field, "drugResistMutations")
				.description("All drug resistance mutations (DRMs) of this drug class."))
			.field(field -> newMutationSet(name, field, "surveilDrugResistMutations")
				.description("All surveillance drug resistance mutations (SDRMs) of this drug class."))
			.field(field -> newMutationSet(name, field, "rxSelectedMutations")
				.description("All treatment selected mutations (TSMs) of this drug class."))
			.field(field -> field
				.type(GraphQLBoolean)
				.name("hasDrugResistMutations")
				.description(
					"Indicate if any mutation is classified as " +
					"Drug Resistance Mutation (DRM) for this drug class."))
			.field(field -> field
				.type(GraphQLBoolean)
				.name("hasSurveilDrugResistMutations")
				.description(
					"Indicate if any mutation is classified as Surveillance " +
					"Drug Resistance Mutation (SDRM) for this drug class."))
			.field(field -> field
				.type(GraphQLBoolean)
				.name("hasRxSelectedMutations")
				.description(
					"Indicate if any mutation is classified as Treatment " +
					"Selected Mutation (TSM) for this drug class."))
			.build()
		)
	);
}
